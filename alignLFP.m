function [LFPmat,Sft] = alignLFP(parentDir,ptID,alignSpot,run1_fName,run2_fName)

% ALIGNLFP aligns LFP data from two files (last two input args)
%   on 'alignSpot': 'stimulus' or 'response'
%
% first output is a tensor of LFP data with dimensions
%       [channels X samples X trials]
%
% second output is a tensor of spectrogram data with dimensions
%       [channels X samples X frequencies X trials]
%

% % !!! delete these if using functionally !!!
% alignSpot = 'response';
% ptID = 'YDX1';
% parentDir = 'D:\Data\Elliot\AlgoPlaceCells\NotBirds_data';
% run1_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-004_task-Birds-Object-run-01_blk-01\EMU-004_subj-YDX_task-Birds-Object_run-01_blk-01_NSP-1.ns3';
% run2_fName = 'D:\Data\Elliot\AlgoPlaceCells\EMU-005_task-Birds-Object-run-01_blk-02\EMU-005_subj-YDX_task-Birds-Object_run-01_blk-02_NSP-1.ns3';

% ptID = 'YEA1';

% loading preprocessed data. Assumes the file name = ptID
load(fullfile(parentDir,[ptID '.mat']))

% then load the LFP.
NSX_1 = openNSx(run1_fName);
NSX_2 = openNSx(run2_fName);
% Concatenating LFP data from both blocks.
tmpData = cat(2,NSX_1.Data{2},NSX_2.Data{2});

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
        [b1,a1] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        tmp(ch,:) = filtfilt(b1,a1,resample(double(tmpData(ch,:)),Fnew,Fs));

        [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
        LFP(ch,:) = filtfilt(b2,a2,tmp(ch,:));
    else
        LFP(ch,:) = resample(double(tmpData(ch,:)),Fnew,Fs);
    end
end
clear NSX_1 NSX_2 NEV tmp tmpData
Fs = Fnew; clear Fnew;

% how much time around each alignment spot you want to include...
pre = 4;
post = 4;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

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

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = nTrials:-1:1
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = LFP(1:nChans,floor(Fs*tTimes(tt))-Fs*pre:floor(Fs*tTimes(tt))+Fs*post);

    % making spectrograms - the code is in this function
    specgram = true;
    if specgram
        for ch2 = nChans:-1:1
            [W,period,scale] = basewaveERP(LFPmat(ch2,:,tt),Fs,1,200,6,0);
            %Sft(ch2,:,:,tt) = abs(W);
            Sft(ch2,:,:,tt) = abs((W))...
                ./repmat(nanmean(nanmean(abs(W(:,tSec>-1 & tSec<-2)),2),3),1,size(W,2),size(W,3));
        end
        % frequencies.
        scaleFreqs = period.^-1;
        fprintf('\ndone generating spectrograms for all channels for trial %d of %d',tt,nTrials)
    end
end

% rule categories
cats = eval([ptID '.stepVarsFlat(:,18);']);
if exist('nanTime','var')
    cats(nanTime) = [];
end

% plotting for categories...





% plotting for the alignment issue...
plt = false;
if plt
    for pl = 1:3
        % plotting a matrix of ERPs.
        subplot(2,3,pl)
        imagesc(tSec,1:nChans,zscore(squeeze(mean(LFPmat(:,:,cats==pl),3)),1,2))
        axis xy

        switch pl
            case 1
                ylabel('channels')
                title('targets')
            case 2
                xlabel('time (s)')
                title('distractors')
            case 3
                title('irrelevants')
        end

        if specgram
            % plotting a grand average spectrogram.
            subplot(2,3,pl+3)
            imagesc(tSec(tSec>-3 & tSec<3),scaleFreqs,squeeze(nanmean(nanmean(Sft(:,:,(tSec>-3 & tSec<3),cats==pl),4))))
            axis xy
            switch pl
            case 1
                ylabel('frequency (Hz)')
                title('targets')
            case 2
                xlabel('time (s)')
                title('distractors')
            case 3
                title('irrelevants')
            end
        else
            % plotting ERPs aligned at zero
            subplot(4,3,pl+6)
            plot(tSec,smoothdata(zscore(squeeze(mean(LFPmat(1:5:end,:,cats==pl),3)),1,2)','movmean',10)'./4)
            axis tight
            xlim([-0.5 1.5])

            % plotting ERPs aligned at -3
            subplot(4,3,pl+9)
            plot(tSec,smoothdata(zscore(squeeze(mean(LFPmat(1:5:end,:,cats==pl),3)),1,2)','movmean',10)'./4)
            axis tight
            xlim([-3.5 -1.5])

            switch pl
            case 1
                ylabel('ERP (uV)')
                title('targets')
            case 2
                xlabel('time (s)')
                title('distractors')
            case 3
                title('irrelevants')
            end
        end
    end
    colormap turbo
    halfMaximize(gcf,'right')

    saveas(gcf,[ptID '_LFPsummary_alignedOn' alignSpot '.pdf'])
end

% save location
saveDir = 'D:\Data\Elliot\AlgoPlaceCells\prePro\';

% for putting in your data structure.
% ... not sure how much detail to include in the details struct.
eval([ptID '.lfp_' alignSpot ' = LFPmat;'])
eval([ptID '.lfp_' alignSpot '_details.Fs = Fs;'])
eval([ptID '.lfp_' alignSpot '_details.preTime = pre;'])
eval([ptID '.lfp_' alignSpot '_details.postTime = post;'])
eval([ptID '.lfp_' alignSpot '_details.tSec = tSec;'])

fprintf('\nsaving data to %s...',saveDir)
if specgram
    eval([ptID '.Sft_' alignSpot ' = Sft;'])
    eval([ptID '.Sft_' alignSpot '_freqs = scaleFreqs;'])
    save([saveDir ptID '_LFPandSft_alignedOn_' alignSpot '.mat'],ptID,'-v7.3')
else
    save([saveDir ptID '_LFPonly_alignedOn_' alignSpot '.mat'],ptID,'-v7.3')
end


end


function [wave,period,scale,coi]=basewaveERP(Y,adrate,loF,hiF,k0,waitc)
% BASEWAVE2 convolves a Morlet-based wavelet with the data to compute the
% time-frequency representation of any signal.  The function returns the
% complex number representation of these signals, and it can be utilized as
% follows:
%
%   [WAVE,PERIOD,SCALE,COI] = BASEWAVE4(DATA,Fs,LOWF,HIF,K0,WAITH)
%
% where the variables are:
%   OUTPUT:
%       WAVE - time-frequency representation of the signal (power and angle can
%          be calculated from these complex numbers);
%       PERIOD - inverse of the frequency scale;
%       SCALE - wavelet decomposition frequencies
%       COI - cone of influence indicates where the wavelet analysis is skewed
%         because of edge artifacts;
%   INPUT:
%       DATA - the signal in the time domain;
%       Fs - sampling frequency;
%       LOWF - lower frequency of the range for which the transform is to be
%          done;
%       HiF - higher frequency of the range for which the transform is to be
%         done;
%       K0 - constant that has to do with converting from scales to fourier
%         space
%       WAITH - a handle to the qiqi waitbar. (usually 0)
%
% See also BASEWAVE, ABRA, WAVELET.

% Original code written by Torrence
% Modified by Peter Lakatos many times
% Modified by Ankoor Shah to return the complex representation
% cleaned up by Elliot Smith (20131126)

dt      = 1/adrate;
dj      = 0.3;                      % adjusts the number of scales. 
s0      = 1/(hiF+(0.1*hiF));    % the smallest resolvable scale
pad     = 0;                        % zero padding (1-on, 0-off)
% k0      = 6;                        % the mother wavelet parameter (wavenumber), default is 6.
frq_lo  = loF;
frq_hi  = hiF;
n = length(Y);
J1 = (log2(n*dt/s0))/dj;                  % [Eqn(10)] (J1 determines the largest scale)

%....waitbar
if waitc==0;
else waitbar(0, waitc);
end

%....az adatsor atalakitasa
n1 = length(Y);
x(1:n1) = Y - mean(Y);
if (pad == 1)
    base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
    x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];                               % k n/2 elembol all
k = k*((2*pi)/(n*dt));
k = [0, k, -k(fix((n-1)/2):-1:1)];              % k ismet n elembol all - angular frequency


%....compute FFT of the (padded) time series (DFT)
f = fft(x);                                     % [Eqn(3)]

%....construct SCALE array
scale = s0*2.^((0:J1)*dj);                      % [Eqn(9)]  (choice of scales to use in wavelet transform ([Eqn(9)]))
fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2));  % Scale-->Fourier [Sec.3h]
coi = fourier_factor/sqrt(2);                   % Cone-of-influence [Sec.3g]
period = fourier_factor*scale;
xxx=min(find(1./period<frq_lo));
period=fliplr(period(1:xxx));
scale=fliplr(scale(1:xxx));
%....find freqency-equivalent scales
% aa=1;
% for a=frq_lo:0.2:frq_hi
%     ii(aa)=max(find(period<(1/a)));
%     scale1(aa)=scale(ii(aa));
%     aa=aa+1;
% end
%
% fscale=aa-1;

scale1=scale;
fscale=size(period,2);

%....construct empty WAVE array
wave = zeros(fscale,n);     % define the wavelet array
wave = wave + 1i*wave;       % make it complex

%....loop through all scales and compute transform
for a1 = 1:fscale
    if waitc==0;
    else waitbar(a1/fscale, waitc);
    end
    expnt = -(scale1(a1).*k - k0).^2/2.*(k > 0.);
    norm = sqrt(scale1(a1)*k(2))*(pi^(-0.25))*sqrt(n);      % total energy=N   [Eqn(7)]
    daughter = norm*exp(expnt);
    daughter = daughter.*(k > 0.);                          % Heaviside step function
    wave(a1,:) = ifft(f.*daughter);                         % wavelet transform[Eqn(4)]
    clear expnt,daughter;
end

period = fourier_factor*scale1;                                 % az ???jfajta, fourier frekvenci???knak megfelelo period
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];   % COI [Sec.3g]

% wave = wave(Laxis:fscale,1:n1);                                  % get rid of padding before returning
% period = period(Laxis:fscale);

%powerx = (abs(wave)).^2;                                        % compute wavelet power spectrum
%powerx = (powerx*1000)/adrate;                                  % ez egy korrekci???, ami a torrence-ben nem volt benne

end


