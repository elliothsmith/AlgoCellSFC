function [SFCdata] = alignLFPandSpikes(parentDir,ptID,alignSpot,run1_fName,run2_fName)
% ALIGNLFPANDSPIKES spike field coherence for algo place cell project.
%   on 'alignSpot': 'stimulus' or 'response'
%

% % !!! delete these if using functionally !!!
% alignSpot = 'response';
% ptID = 'YDX1';
% parentDir = 'D:\Data\Elliot\AlgoPlaceCells\NotBirds_data';
% run1_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-004_task-Birds-Object-run-01_blk-01\EMU-004_subj-YDX_task-Birds-Object_run-01_blk-01_NSP-1.ns3';
% run2_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-005_task-Birds-Object-run-01_blk-02\EMU-005_subj-YDX_task-Birds-Object_run-01_blk-02_NSP-1.ns3';

% ptID = 'YEA1';

% loading preprocessed data and counting units.
% Assumes the file name = ptID
load(fullfile(parentDir,[ptID '.mat']))
eval(['unitData = ' ptID '.psths;'])
nUnits = length(unitData);

% then load the LFP.
NSX_1 = openNSx(run1_fName);
NSX_2 = openNSx(run2_fName);
% Concatenating LFP data from both blocks.
tmpData = cat(2,NSX_1.Data{2},NSX_2.Data{2});

% macroelectrode labels. (woof!)
macroLabels = {NSX_1.ElectrodesInfo.Label}';
microLabels = {unitData.electrodeLabel}';

% data parameters
nSamps = size(tmpData,2);
nChans = size(tmpData,1);
% sampling rates
Fs = NSX_1.MetaTags.SamplingFreq;

% resampling LFP at Fnew sampling frequency and notch filtering...
notchFilter = true;
Fnew = 400;
for ch = nChans:-1:1
    if notchFilter
        [b,a] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        tmp(ch,:) = filtfilt(b,a,resample(double(tmpData(ch,:)),Fnew,Fs));

        [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
        LFP(ch,:) = filtfilt(b2,a2,tmp(ch,:));
    else
        LFP(ch,:) = resample(double(tmpData(ch,:)),Fnew,Fs);
    end
end
clear NSX_1 NSX_2 NEV tmp tmpData
Fs = Fnew; clear Fnew;

% how much time around each alignment spot you want to include...
pre = 3;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);
SFCdata.pre = pre;
SFCdata.post = post;

% timing
if isequal(alignSpot,'stimulus')
    tTimes = eval([ptID '.stepVarsFlat(:,7);']);
elseif isequal(alignSpot,'response')
    tTimes = eval([ptID '.stepVarsFlat(:,11);']);
    % one of the response times is a NaN???
    nanTime = isnan(tTimes);
    tTimes(nanTime) = [];
end
nTrials = length(tTimes);

% [20220417] I originally wrote this code to look at trial by trial
% comparisons, but I am starting to think that doesn't make much sense, so
% just going to do trial averaged stuff here.

% first figure out which categories to look at...
comparison = 'targetstatus'; % 'stimulusidentity' 'targetstimulus'
switch comparison
    case {'stimulusidentity'}
        % stimulus identity==6
        cats = eval([ptID '.stepVarsFlat(:,6);']);
    case {'targetstatus'}
        % traget status ==18
        cats = eval([ptID '.stepVarsFlat(:,18);']);
    case {'targetstimulus'}
        % target stimulus ==23
        cats = eval([ptID '.stepVarsFlat(:,23);']);
end
% saving
SFCdata.comparison = comparison; 

if exist('nanTime','var')
    cats(nanTime) = [];
end

% looping over trials to epoch data.
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = nTrials:-1:1
    % epoch the LFP data here [channels X samples X trials]
    LFPmat(:,:,tt) = LFP(1:nChans,floor(Fs*tTimes(tt))-Fs*pre:floor(Fs*tTimes(tt))+Fs*post);

    % looping over units (and trials).
    for un= 1:nUnits
        % put spikes in structure
        spikes.unit(un).trial(tt).times = unitData(un).spikeTimes(unitData(un).spikeTimes>tTimes(tt)-pre & unitData(un).spikeTimes<tTimes(tt)+post)... - (tTimes(tt)-preZ)
            - repmat(tTimes(tt)-pre,length(unitData(un).spikeTimes(unitData(un).spikeTimes>tTimes(tt)-pre & unitData(un).spikeTimes<tTimes(tt)+post)),1);
    end
end
% now I have the mfing data in chroniux format, I can generate
% cohereograms for each trial type and do permutation tests.

% first have to specify chronux params.
params.tapers = [9 17];
params.pad = 0;
params.Fs = Fs;
params.fpass = [1 50];
params.err = [1 0.05];
params.trialave = 1;
movingWin = [1 0.05];
% not the possibility of doing two different amounts of smoothing for
% different frequecny bands, as in Womelsdorf et al. (2006).

% saving params
SFCdata.params = params;

% starting parallel pool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
end
pul = parpool(16);

% looping over units to examine coherence between each units FR and the LFP
% on the nearest macroContact.
for un = 1:nUnits

    % figuring out which macro to look at.
    % starting by constructing the macro label from the micro label.
    nearestMacroIdx = contains(macroLabels,[deblank(microLabels{un}(2:end-2)) '01']);
    macroChan = find(nearestMacroIdx);
    nearestMacroLabel = macroLabels{nearestMacroIdx};

    % saving
    SFCdata.unit(un).nearestMacroLabel = nearestMacroLabel;

    nCats = length(unique(cats));
    figure(un)
    for ct = 1:nCats
        % calculate coherence
        [C,~,~,~,~,t,f,~,~,~]=cohgramcpt(squeeze(LFPmat(macroChan,:,cats==ct)),spikes.unit(un).trial(cats==ct),movingWin,params,0);

        % do permutation testing between categories. PARFOR
        tic
        nPerms = 200;
        parfor prm = 1:nPerms
            spikePerms = randperm(length(cats),sum(cats==ct));
            LFPPerms = randperm(length(cats),sum(cats==ct));
            % regenerate cohereograms for random trial permutations
            [permC,~,~,~,~,~,~,~,~,~]=cohgramcpt(squeeze(LFPmat(macroChan,:,LFPPerms)),spikes.unit(un).trial(spikePerms),movingWin,params,0);

            % counting permutation values less than C
            gtC(prm,:,:) = C>permC;

            % retaining cluster sizes
            clusterSizesPRM{prm} = cell2mat(cellfun(@length,bwboundaries(squeeze(gtC(prm,:,:))),'UniformOutput',false));
        end
        
        permClusterSizes = [];
        for zz = 1:nPerms
            permClusterSizes = cat(1,permClusterSizes,clusterSizesPRM{zz});
        end
        criticalClusterSize = min(maxk(permClusterSizes,floor(0.05*length(permClusterSizes))));

        % finding significant clusters and their sizes
        criticalNum = floor(0.95*nPerms);
        sigMask = squeeze(sum(gtC))>criticalNum;
        cc = bwconncomp(sigMask);
        bb = bwboundaries(sigMask);
%         sigClusterInfo = regionprops(cc,'Area','BoundingBox'); 

        % cluster correction
        clusterSizes = cell2mat(cellfun(@length,bb,'UniformOutput',false));
        sigClusters = clusterSizes>criticalClusterSize;    
        sigClusterIdcs = find(sigClusters);
        nSigClusters = sum(sigClusters);

        % time relative to alingment.
        tSec = t-pre;

        % reporting elapsed time.
        A = toc;
        fprintf('\n\npermutation test with %d permutations took %.2f minutes to run...\n',nPerms,A./60)

        % plot the results...
        subplot(nCats,1,ct)
        hold on

        % plot coherogram
        imagesc(tSec,f,C')
        line([0 0],[1 50],'linestyle','--','color','k')
        % plot bounding box around significant clusters
        boundaryTF = {};
        for cls = 1:nSigClusters
%             BB = sigClusterInfo(sigClusterIdcs(cls)).BoundingBox;
%             rectangle('Position',[tSec(ceil(BB(1))),f(ceil(BB(2)))-(f(BB(4))./2),t(BB(3)),f(BB(4))./2],'EdgeColor','k','LineWidth',2);
              boundary = bb{sigClusterIdcs(cls)};
              plot(tSec(boundary(:,1)),f(boundary(:,2)),'k-.')

              % you may want to save the boundaries of significant
              % clusters. 
              boundaryTF{cls} = [tSec(boundary(:,1)) f(boundary(:,2))];
        end

        hold off
        xlabel('time (s)')
        ylabel('frequecny (Hz)')
        title([comparison ': ' num2str(ct)])
        axis square xy tight
        colorbar

        % saving cohgrams
        SFCdata.unit(un).time = tSec;
        SFCdata.unit(un).freq = f;
        SFCdata.unit(un).category(ct).cohereogram = C;
        SFCdata.unit(un).category(ct).numSigClusters = nSigClusters;
        SFCdata.unit(un).category(ct).clusterTimeFreqBoundary = boundaryTF;


    end

    % saving figure
    halfMaximize(un,'left')
    saveas(gcf,['D:\Figs\Elliot\AlgoPlaceCells\' ptID '_unit' num2str(un) '_' microLabels{un} '_' comparison '_coherograms_alignedOn' alignSpot '.pdf'])
    close(un)

end

% saving data. 
save(['D:\Data\preProcessed\algoPlaceCells\' ptID '_' microLabels{un} '_' comparison '_coherograms_alignedOn' alignSpot '.pdf'],'SFCdata','-v7.3')














%     % making spectrograms [channels X time X freq X trials] - the code is in this function
%     for ch2 = nChans:-1:1
%         [W,period,~] = basewaveERP(LFPmat(ch2,:,tt),Fs,1,200,6,0);
%         %Sft(ch2,:,:,tt) = abs(W);
%         Sft(ch2,:,:,tt) = abs((W))...
%             ./repmat(nanmean(nanmean(abs(W(:,tSec>-1 & tSec<-2)),2),3),1,size(W,2),size(W,3));
%     end
%     % frequencies.
%     scaleFreqs = period.^-1;
%     fprintf('\ndone generating spectrograms for all channels for trial %d of %d',tt,nTrials)
%
%
%
%
%     % Now generating psths for each scale to do coherence.
%     nFs = length(scaleFreqs);
%     nUnits = length(unitData);
%     % looping over units
%     for un= 1:nUnits
%         % put spikes in structure
%         spikes.unit(un).trial(tt).times = unitData(un).spikeTimes(unitData(un).spikeTimes>tTimes(tt)-pre & unitData(un).spikeTimes<tTimes(tt)+post) - (tTimes(tt)-preZ)
%            % - repmat(tTimes(tt)-pre,length(unitData(un).spikeTimes(unitData(un).spikeTimes>tTimes(tt)-pre & unitData(un).spikeTimes<tTimes(tt)+post)),1);
%
%         % psths with differen kernel widths for SFC calculation.
%         for ff = 1:nFs
%             kernelWidth = period(ff)./4;
%
%             % calculating averaged psth to get time base...
%             [~,~,t] = psthBins(spikes.unit(un).trial(1).times, kernelWidth,1,1,pre+post);
%
%             % single trial spike periodograms [unit X frequency X trials X time]
%             if isempty(spikes.unit(un).trial(tt).times)
%                 Rft(un,ff,tt,:) = zeros(1,1,1,length(t));
%             else
%                 [~,Rft(un,ff,tt,:),t] = psthBins(spikes.unit(un).trial(tt).times, kernelWidth,1,1,pre+post);
%             end
%
%             % timing
%             tSecCue = linspace(-pre,post,length(t));
%         end
%     end
% end



% keyboard
%
% % rule categories
% cats = eval([ptID '.stepVarsFlat(:,18);']);
% if exist('nanTime','var')
%     cats(nanTime) = [];
% end
%
% % plotting for categories...
%
%
% % plotting for the alignment issue...
% plt = false;
% if plt
%     for pl = 1:3
%         % plotting a matrix of ERPs.
%         subplot(2,3,pl)
%         imagesc(tSec,1:nChans,zscore(squeeze(mean(LFPmat(:,:,cats==pl),3)),1,2))
%         axis xy
%
%         switch pl
%             case 1
%                 ylabel('channels')
%                 title('targets')
%             case 2
%                 xlabel('time (s)')
%                 title('distractors')
%             case 3
%                 title('irrelevants')
%         end
%
%         if SFC
%             % plotting a grand average spectrogram.
%             subplot(2,3,pl+3)
%             imagesc(tSec(tSec>-3 & tSec<3),scaleFreqs,squeeze(nanmean(nanmean(Sft(:,:,(tSec>-3 & tSec<3),cats==pl),4))))
%             axis xy
%             switch pl
%                 case 1
%                     ylabel('frequency (Hz)')
%                     title('targets')
%                 case 2
%                     xlabel('time (s)')
%                     title('distractors')
%                 case 3
%                     title('irrelevants')
%             end
%         else
%             % plotting ERPs aligned at zero
%             subplot(4,3,pl+6)
%             plot(tSec,smoothdata(zscore(squeeze(mean(LFPmat(1:5:end,:,cats==pl),3)),1,2)','movmean',10)'./4)
%             axis tight
%             xlim([-0.5 1.5])
%
%             % plotting ERPs aligned at -3
%             subplot(4,3,pl+9)
%             plot(tSec,smoothdata(zscore(squeeze(mean(LFPmat(1:5:end,:,cats==pl),3)),1,2)','movmean',10)'./4)
%             axis tight
%             xlim([-3.5 -1.5])
%
%             switch pl
%                 case 1
%                     ylabel('ERP (uV)')
%                     title('targets')
%                 case 2
%                     xlabel('time (s)')
%                     title('distractors')
%                 case 3
%                     title('irrelevants')
%             end
%         end
%     end
%     colormap turbo
%     halfMaximize(gcf,'right')
%
%     saveas(gcf,[ptID '_LFPsummary_alignedOn' alignSpot '.pdf'])
% end
%
% % save location
% saveDir = 'D:\Data\Elliot\AlgoPlaceCells\prePro\';
%
% % for putting in your data structure.
% % ... not sure how much detail to include in the details struct.
% eval([ptID '.lfp_' alignSpot ' = LFPmat;'])
% eval([ptID '.lfp_' alignSpot '_details.Fs = Fs;'])
% eval([ptID '.lfp_' alignSpot '_details.preTime = pre;'])
% eval([ptID '.lfp_' alignSpot '_details.postTime = post;'])
% eval([ptID '.lfp_' alignSpot '_details.tSec = tSec;'])
%
% fprintf('\nsaving data to %s...',saveDir)
% if SFC
%     eval([ptID '.Sft_' alignSpot ' = Sft;'])
%     eval([ptID '.Sft_' alignSpot '_freqs = scaleFreqs;'])
%     save([saveDir ptID '_LFPandSft_alignedOn_' alignSpot '.mat'],ptID,'-v7.3')
% else
%     save([saveDir ptID '_LFPonly_alignedOn_' alignSpot '.mat'],ptID,'-v7.3')
% end
%
%
% end
