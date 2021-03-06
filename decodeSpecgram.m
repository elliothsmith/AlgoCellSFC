% load data
ptID = 'YEA1';
alignSpot = 'stimulus';
saveDir = 'D:\Data\Elliot\AlgoPlaceCells\prePro\';
load([saveDir ptID '_LFPandSft_alignedOn_' alignSpot '.mat'],ptID,'-v7.3')

% vectorize spectrograms
eval(['Sft = ' ptID '.Sft_' alignSpot ';'])
eval(['f = ' ptID '.Sft_' alignSpot '_freqs;'])


% start parpool

% run habiba's code:

% options.
taskVar = 'targetStimulus'; 
groupVar = '';  % a grouping variable, to decode within categories of this secondary variable; takes all the same options as taskVar
numGroups = 3;  % resets to 1 if no grouping variable is used
shuffleMethod = 'sessionPerm';  % 'sessionPerm'(permutes session information), 'behaviorPerm' (randomly permutes behavior labels) and 'circular' (circularly shuffles spike times)
combinedSessionsFlag = 0;   % whether to combine across all behavior sessions to create a pseudo-population; doesn't make sense to set this to 1 when sessionNum has been specified; maybe in future I could expand this option to combine across only specific sessions.
decoderType = 'lda';    % 'lda', 'logistic', 'svm', 'knn'
params = struct('gamma', 0);  % for LDA; not currently used otherwise
% params = [];
alignToEvent = 'stimulus';  % 'stimulus' or 'response'
timeWindow = {[-1 -0.5], [-0.5 0], [0 0.5]}; % s; used as separate features for each neuron
% timeWindow = {[-1 2]};    % I use this when slidingWindowFlag = true
% timeWindow = {[0.2 0.7]}; % I use this if I only want to use one time feature
slidingWindowFlag = false;  % if true, decodes on a sliding time window, specified above
slidingWindowSize = 0.5;        % s; only relevant if slidingWindowFlag = true
slidingWindowIncrement = 0.25;   % s; only relevant if slidingWindowFlag = true

%% figure number index; increments as new figures are created
figNum = 0;

%% decode on real data
disp('Decoding on real data');
fraction_correct = [];  % # of sessions x # of groups x # of folds
confusion_matrix = [];  % # of sessions x # of groups x # of categories x # of categories
all_models = {};    % # of sessions x # of groups x # of folds
classLabels = {};   % # of sessions x # of groups (in case these are different for different sessions, which mostly they shouldn't be)
tic
for s = sessionIxes
    % select appropriate data for session (or all sessions combined)
    if combinedSessionsFlag
        if slidingWindowFlag
            X = S.X_full;   % this includes full PSTH
        else
            X = S.X;        % this includes only averaged firing rates over the desired time window
        end
        Y = S.Y;
        G = S.G;
    else
        if slidingWindowFlag
            X = S(s).X_full;
        else
            X = S(s).X;
        end
        Y = S(s).Y;
        G = S(s).G;
    end
    for gix = 1:numGroups
        Xg = X(G==gix,:,:);     % this line will work even if X is only 2D
        Yg = Y(G==gix);
        classLabels{s, gix} = unique(Yg);
        if slidingWindowFlag    % will iterate through time windows and decode
            timeWindow_flat = timeWindow{1};
            timeBins_data = round(S.timeBins_s{1}, 2);
            timeBins_start = round(timeWindow_flat(1) : slidingWindowIncrement : timeWindow_flat(end) - slidingWindowSize, 2);
            timeBins_end = round(timeWindow_flat(1) + slidingWindowSize : slidingWindowIncrement : timeWindow_flat(end), 2);
            parfor (tix = 1:length(timeBins_start), numWorkers) 
        %         disp([num2str(tix) '/' num2str(length(timeBins_start))]);
                tix_start = find(timeBins_data == timeBins_start(tix));
                tix_end = find(timeBins_data == timeBins_end(tix));
                [fraction_correct(s,gix,tix,:), confusion_matrix(s,gix,tix,:,:)] = decodeData(Xg(:,:, tix_start:tix_end), Yg, decoderType, params);
        %         disp([num2str(tix) '/' num2str(length(timeBins_start)) ' done']);
            end
        else
            [fraction_correct(s,gix,:), confusion_matrix(s,gix,:,:), all_models(s,gix,:)] = decodeData(Xg, Yg, decoderType, params);
        end
    end
end
toc


%% visualize average weights from classifiers across all folds
% % (this bit of code only works in a subset of cases; haven't fleshed it out
% % for the rest yet)
% if combinedSessionsFlag && ~slidingWindowFlag && numGroups == 0
%     figNum = figNum + 1; f = figure(figNum); set(f, 'Position', [50 50 1000 500]); clf; hold on;
%     numClasses = size(confusion_matrix,3);
%     numNeurons = size(X,3);
%     numFolds = size(fraction_correct,3);
%     
%     if strcmpi(decoderType, 'LDA')
%         classifierWeights = NaN(numClasses, numClasses, numNeurons, numFolds);
% 
%         for c1 = 1:numClasses
%             for c2 = (c1+1):numClasses
%                 for cv = 1:numFolds
%                     if c1 == c2
%                         continue;
%                     end
%                     classifierWeights(c1, c2, :, cv) = all_models{1,cv}.Coeffs(c1,c2).Linear;
%                 end
%                 bar(1:numNeurons, squeeze(median(classifierWeights(c1, c2, :, :), 4)), 'displayname', [num2str(c1) '/' num2str(c2)]);
%             end
%         end
%     elseif strcmpi(decoderType, 'logistic')
%         classifierWeights = NaN(numClasses, numNeurons, numFolds);
%         for c1 = 1:numClasses
%             for cv = 1:numFolds
%                 classifierWeights(c1,:,cv) = all_models{1,cv}.BinaryLearners{c1}.Beta;
%             end
%             bar(1:numNeurons, squeeze(median(classifierWeights(c1,:,:),3)), 'displayname', num2str(c1));
%         end
%     end
%     legend('show');
%     xlabel('Neuron #');
%     ylabel('Classifier weight');
%     saveas(gca, [outputDir decoderType '_weights.png']);
%     clearvars f;
% end

% %% save all variables to output folder under workspace.mat
% save([outputDir 'workspace.mat'], ...
%     'alignToEvent', 'all_models', 'classLabels', 'combinedSessionsFlag', 'sessionNum', ...
%     'confusion_matrix', 'dataFile', 'dataPath', 'decoderType', 'fraction_correct', ...
%     'groupVar', 'name', 'numGroups', 'numPerms', 'numSessions', 'outputDir', ...
%     'params', 'region', 'shuffleMethod', 'slidingWindowFlag', 'slidingWindowIncrement', 'slidingWindowSize', ...
%     'taskVar', 'timeStr', 'timeWindow', '-v7.3');

%% do the same for shuffled data
disp('Decoding on shuffled data');
fraction_correct_shuff = [];    % # of sessions x # of groups x # of permutations (x # of time bins, if sliding window) x # of folds
confusion_matrix_shuff = [];    % # of sessions x # of groups x # of permutations (x # of time bins, if sliding window) x # of categories x # of categories
% all_models_shuff = {};    % takes up a ton of space; use at your own peril!

for p = 1:numPerms
    tic
    % display purposes
    if mod(p,numPerms/10)==0
        disp(['*' num2str(p) '/' num2str(numPerms)]);
    end
    if combinedSessionsFlag
        if ~isempty(groupVar)
            permDataFile = ['./formatDataForDecoding_' taskVar 'X' groupVar '_' region '_' alignToEvent '_time' timeStr '/pseudopop_collapsed_' shuffleMethod num2str(p) '.mat'];
        else
            permDataFile = ['./formatDataForDecoding_' taskVar '_' region '_' alignToEvent '_time' timeStr '/pseudopop_collapsed_' shuffleMethod num2str(p) '.mat'];
        end
        if ~forceRecompute && exist(permDataFile, 'file')
%             disp('Loading permuted data');
            S_shuff = load(permDataFile); S_shuff = S_shuff.cData;
        else
            [S_shuff,~] = formatDataForDecoding(dataPath, taskVar, groupVar, region, alignToEvent, timeWindow, shuffleMethod, 1, num2str(p));
        end
    else
        if ~isempty(groupVar)
            permDataFile = ['./formatDataForDecoding_' taskVar 'X' groupVar '_' region '_' alignToEvent '_time' timeStr '/pseudopop_' shuffleMethod num2str(p) '.mat'];
        else
            permDataFile = ['./formatDataForDecoding_' taskVar '_' region '_' alignToEvent '_time' timeStr '/pseudopop_' shuffleMethod num2str(p) '.mat'];
        end
        if ~forceRecompute && exist(permDataFile, 'file')
            S_shuff = load(permDataFile); S_shuff = S_shuff.sData;
        else
            [~, S_shuff] = formatDataForDecoding(dataPath, taskVar, '', region, alignToEvent, timeWindow, shuffleMethod, 1, num2str(p));
        end
    end
    
    for s = sessionIxes
        if combinedSessionsFlag
            if slidingWindowFlag
                X_shuff = S_shuff.X_full;
            else
                X_shuff = S_shuff.X;
            end
            Y_shuff = S_shuff.Y;
            G_shuff = S_shuff.G;
        else
            if slidingWindowFlag
                X_shuff = S_shuff(s).X_full;
            else
                X_shuff = S_shuff(s).X;
            end
            Y_shuff = S_shuff(s).Y;
            G_shuff = S_shuff(s).G;
        end
        for gix = 1:numGroups
            Xg_shuff = X_shuff(G_shuff==gix,:,:);         % this line will work even if X is only 2D
            Yg_shuff = Y_shuff(G_shuff==gix);
            if slidingWindowFlag
                temp_fraction_correct_shuff = [];   % I had a parfor somewhere, which is why I used these, but can't remember which loop it was now... 
                temp_confusion_matrix_shuff = [];
                parfor (tix = 1:length(timeBins_start), numWorkers)
        %             disp([num2str(tix) '/' num2str(length(timeBins_start))]);
                    tix_start = find(timeBins_data == timeBins_start(tix));
                    tix_end = find(timeBins_data == timeBins_end(tix));
                    if isempty(tix_start) || isempty(tix_end)
                        keyboard;
                    end
                    [temp_fraction_correct_shuff(tix,:), temp_confusion_matrix_shuff(tix,:,:), ~] = decodeData(Xg_shuff(:,:, tix_start:tix_end), Yg_shuff, decoderType, params);
                end
                fraction_correct_shuff(s,gix,p,:,:) = temp_fraction_correct_shuff;
                confusion_matrix_shuff(s,gix,p,:,:,:) = temp_confusion_matrix_shuff;
%                 all_models_shuff(s,gix,p,:) = temp_all_models_shuff;
            else
                [fraction_correct_shuff(s,gix,p,:), confusion_matrix_shuff(s,gix,p,:,:), ~] = decodeData(Xg_shuff, Yg_shuff, decoderType, params);
            end
        end
    end
    if mod(p, numPerms/10)==0
        try
            %% incrementally save results, in case everything crashes... 
            save([outputDir 'workspace.mat'], ...
                'confusion_matrix_shuff', 'fraction_correct_shuff', '-append');
        catch
            fclose('all');
        end
    end
    toc;
end
%% save results at the end
save([outputDir 'workspace.mat'], ...
    'confusion_matrix_shuff', 'fraction_correct_shuff', '-append');

%% visualize results
for s = sessionIxes
    if ~combinedSessionsFlag
        sessionString = ['s' num2str(s) '_'];
    else
        sessionString = '';
    end
    if slidingWindowFlag
        figNum = figNum + 1; figure(figNum); clf; 
        for gix = 1:numGroups
            subplot(numGroups, 1, gix); hold on;
            % plot average cross-validated decoding accuracy on data across time
            plot(timeBins_start + (slidingWindowSize/2), squeeze(mean(fraction_correct(s,gix,:,:),4)));
            % mean accuracy (over all permutations and all folds) from shuffled data
            plot(timeBins_start + (slidingWindowSize/2), squeeze(mean(mean(fraction_correct_shuff(s,gix,:,:,:),5),3)), 'k');
            % 5th and 95th percentile accuracy on shuffled data (through time)
            plot(timeBins_start + (slidingWindowSize/2), squeeze(prctile(mean(fraction_correct_shuff(s,gix,:,:,:),5),5,3)), '--k');
            plot(timeBins_start + (slidingWindowSize/2), squeeze(prctile(mean(fraction_correct_shuff(s,gix,:,:,:),5),95,3)), '--k');
            % plot max-corrected average cross-validated decoding accuracy on shuffled
            % data across time
            maxTestThreshold = prctile(squeeze(max(mean(fraction_correct_shuff(s,gix,:,:,:),5),[],4)),95);
            hline(maxTestThreshold, 'Color', 'r');
            vline(0);
            xlabel(['Time from ' alignToEvent ' onset (s)']);
            ylabel('Decoding accuracy');
            if numGroups > 1
                title(['Group ' num2str(gix)]);
            end
        end
        saveas(gca, [outputDir sessionString 'accuracyThroughTime.png']);
    else
        figNum = figNum + 1; f = figure(figNum); clf;
        if numGroups > 1
            set(f, 'Position', [50 50 1000 500]);
        end
        for gix = 1:numGroups
            subplot(1, numGroups, gix); hold on;
            % histogram of decoder accuracy on shuffled data
            histogram(squeeze(mean(fraction_correct_shuff(s,gix,:,:),4)));
            % decoder accuracy on real data
            vline(mean(fraction_correct(s,gix,:)));
            % p-value in title
            if numGroups > 1
                title(['Group ' num2str(gix) ': p = ' num2str(mean(squeeze(mean(fraction_correct_shuff(s,gix,:,:),4)) > mean(fraction_correct(s,gix,:))))]);
            else
                title(['p = ' num2str(mean(squeeze(mean(fraction_correct_shuff(s,gix,:,:),4)) > mean(fraction_correct(s,gix,:))))]);
            end
            xlabel('Decoder accuracy');
        end
        saveas(gca, [outputDir sessionString 'accuracy.png']);
        
        figNum = figNum + 1; f = figure(figNum); clf;
        if numGroups > 1
            set(f, 'Position', [50 50 1500 500]);
        end
        for gix = 1:numGroups
            subplot(1, numGroups, gix);
            imagesc(squeeze(confusion_matrix(s,gix,:,:))); colorbar;
            set(gca, 'XTick', classLabels{s,gix}', 'YTick', classLabels{s,gix}')
            title('real data', 'FontSize', 16)
            xlabel('true label', 'FontSize', 16)
            ylabel('predicted label', 'FontSize', 16)
            title(['Group ' num2str(gix)]);
        end
        saveas(gca, [outputDir sessionString 'confusionMatrix.png']);
    end
end
