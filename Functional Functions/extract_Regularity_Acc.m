function [WB_SDMO] = extract_Regularity_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_Regularity_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Regularity measures based on Autocorrelation signal - ACCELEROMETRY (VT and AP)
% * Ref 1) % According to: "Moe-Nilssen, Rolf, and Jorunn L. Helbostad.
% "Estimation of gait cycle characteristics by trunk accelerometry."
% Journal of biomechanics 37, no. 1 (2004): 121-126."
% * Ref 2) Authors: Rispens SM, van Schooten KS, Pijnappels M, Daffertshofer A, Beek PJ, van Dieën JH.
% "Identification of Fall Risk Predictors in Daily Life Measurements"
% Neurorehabil Neural Repair. 2015 Jan;29(1):54-61
%
%% Input
% 1) imu: contains 3D acc and 3D angular velocities (if available) with the following fields:
% imu.acc_AP
% imu.acc_ML
% imu.acc_V
% imu.gyro_roll
% imu.gyro_pitch
% imu.gyro_yaw
% 2) fs: sample rate
% 3) WB_WBD: walking bout detection info (output from WBDetection function with the following fields:
% WB_WBD.start
% WB_WBD.end
% WB_WBD.stepsNumber (not included in this example)
% WB_WBD.TempEvents (not included in this example)
% WB_WBD.breaks (not included in this example)
% WB_WBD.turnings
% WB_WBD.turnings.start
% WB_WBD.turnings.end
% WB_WBD.straightWalk
% WB_WBD.straightWalk.start
% WB_WBD.straightWalk.end
% WB_WBD.straightWalk.stepsNumber (we chose 20 steps in each straight walk just for the example)
% 4) WB_SDMO: Secondary Digital Mobility Outcomes structure (with all features)
% Remark: Initialized as an empty structure
% *** If pre-setted:
% 5) Settings for Spatiotemporal extraction
%       Settings.MinNumEvents = 3;                  %CONSENSUS!             % minimal number of steps to consider the step-to-stpe gait features
%       Settings.GravityConstant = 9.81; % [m/s2]
%       Settings.Fc_LowPassFilter = 20*fs/100;                              % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.Fc_HighPassFilter = 0.1*fs/100;                            % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.Fc_HighPassFilter_AutoCorr = 0.8*fs/100;                   % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.StepDetectMcCameley_EMA = 1;                               % use McCamely algorithm with Encarna modifications (1) or with BAM original (0)
%       Settings for frequency analysis
%       Settings.N_Harm = 20;                                               % number of harmonics used for harmonic ratio, index of harmonicity and phase fluctuation
%       Settings.MaxRangeFrequencyPSD = 10*fs/100;                          % maximal range of frequencies to find frequency of the peak power spectral density (PSD)
%       Settings.MinRangeFrequencyPSD = 0.3*fs/100;                         % miniaml range of frequencies to find frequency of the peak power spectral density (PSD)
%       Settings.StrideDurationRange = [0.4 4.0];                           % Range to search for stride duration (seconds)
%       PlotSettings.SpectralFatures = 1;                                   % activate plots to help understangind the calculation of power spectral density of signals

%% Output
% WB_SDMO(WB_index).StepRegularity_Acc.(NameChannels{iChannel}) = regularity of steps based on autocorrelation signal;
% WB_SDMO(WB_index).StrideRegularity_Acc.(NameChannels{iChannel}) = regularity of strides based on autocorrelation signal;
% WB_SDMO(WB_index).Symmetry_AutocorrelationRatio.(NameChannels{iChannel}) = asymmetry of step-stride (ratio) based on autocorrelation signal;
% WB_SDMO(WB_index). Symmetry_AutocorrelationDiff.(NameChannels{iChannel}) = asymmetry of step-stride (difference) based on autocorrelation signal;

%% Remarks
% *** Comments located on the right side, based on Vrije Universiteit algorithm (original author: Sietse Rispens)
% *** Author: Encarna Micó Amigo, Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
%             Chris Buckley, Contact: Christopher.Buckley2@newcastle.ac.uk

%% History
% 2019/6th/September functionized - Encarna Micó Amigo


%% SETTINGS

if nargin < 5                                                               % if the settings are not predefined:
    % Settings for Spatiotemporal extraction
    Settings.MinNumEvents = 3;                  %CONSENSUS!                 % minimal number of steps to consider the step-to-stpe gait features
    Settings.GravityConstant = 9.81; % [m/s2]
    Settings.Fc_LowPassFilter = 20*fs/100;                                  % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter = 0.1*fs/100;                                % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter_AutoCorr = 0.8*fs/100;                       % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.StepDetectMcCameley_EMA = 1;                                   % use McCamely algorithm with Encarna modifications (1) or with BAM original (0)
    % Settings for frequency analysis
    Settings.N_Harm = 20;                                                   % number of harmonics used for harmonic ratio, index of harmonicity and phase fluctuation
    Settings.MaxRangeFrequencyPSD = 10*fs/100;                              % maximal range of frequencies to find frequency of the peak power spectral density (PSD)
    Settings.MinRangeFrequencyPSD = 0.3*fs/100;                             % miniaml range of frequencies to find frequency of the peak power spectral density (PSD)
    Settings.StrideDurationRange = [0.4 4.0];                               % Range to search for stride duration (seconds)
    PlotSettings.SpectralFatures = 1;                                       % activate plots to help understangind the calculation of power spectral density of signals
end

%% LOOP OVER WALKING BOUTS

%WB_Number = size(WB_WBD,2);                                               %number of walking bouts for loop
%for WB_index = 1:WB_Number                                                %loop over each WB
    
% Info about straight line walks within WB
%StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);              %number/amount of straight walk episodes
straightWalk_start = cell2mat({WB_WBD.start});                             %start of straight WB
straightWalk_end = cell2mat({WB_WBD.end});                                 %end of straight WB
        
% Initialization
StepRegularity = [];
StrideRegularity = [];
Symmetry_AutocorrelationRatio = [];
Symmetry_AutocorrelationDiff = [];
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB

% Definition of signal: start-end
Acc = zeros(size(imu.acc_V(straightWalk_start:straightWalk_end),1),4);
Acc(:,1) = imu.acc_V(straightWalk_start:straightWalk_end);
Acc(:,2) = imu.acc_ML(straightWalk_start:straightWalk_end);
Acc(:,3) = imu.acc_AP(straightWalk_start:straightWalk_end);
if mean(Acc(:,1)) < 0                                                      % signal negative = VT acceleration input is downwards
    Acc(:,1) = -Acc(:,1);                                                  % we invert the signal (to get an upwards direction). If it would be directly provided as upwards, there would not be need to invert it
    Acc(:,3) = -Acc(:,3);                                                  % we invert the signal (to get an forwards direction). If it would be directly provided as forwards, there would not be need to invert it
end %although there is no need no invert it

% Detrend & filter signals
Acc_Detr = detrend(Acc);                                                   % remove the trend/offset of gravity
Acc_Detr_LPFilt = WintFilt_low(Acc_Detr-mean(Acc_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal
Acc_Detr_LPFilt(:,4) = sqrt(sum(Acc_Detr_LPFilt(:,1:3)'.^2)');             %resultant or combined signal

% Initial temporal features calculation
%if ~isfield(WB_WBD,'TempEvents')                                          %not provided as an input (FAU & EPFL work)
AccVT = Acc(:,1);
if Settings.StepDetectMcCameley_EMA
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
else
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
end
IC = TemporalFeatures.InitialContactSample;
StepsNumber = size(IC,2);
% else                                                                % provided as an input
%     IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
%     StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);
% end

% Calculate Regularity based on Autocorrelation signal - Acc
NameChannels = {'VT','ML','AP','Combined'};

% Calculate Power Spectral Density - Acc
if StepsNumber > Settings.MinNumEvents
    % Stride time and regularity from auto correlation
    [Autocorr4x4,Lags] = xcov(Acc_Detr_LPFilt,'unbiased');          % auto-covariance of Moe-Nilssen. Consider also to use xcorr instead: The functions xcorr and xcov estimate the cross-correlation and cross-covariance sequences of random processes. The xcov function estimates autocovariance and cross-covariance sequences. This function has the same options and evaluates the same sum as xcorr, but first removes the means of x and y. "Unbiased": non-attenuated signal
    Autocorr4x4 = Autocorr4x4(Lags >= 0,[1 6 11 16]);               % we get positive lags signal (Lags>=0), we select channels 1, 6, 11 and 16 because the function has as an output the correlation of each of the signals with each of the signals, thus, to get the "auto"correlation, we select particular axes of the final matrix
    AutocorrSum = sum(Autocorr4x4(:,1:3),2);                        % this sum is independent of sensor re-orientation, as long as axes are kept orthogonal. Notice that this is different from the autocorrelation of the resultant/combined signal (4th channel of Autocorr4x4)
    AutocorrAll_ = [Autocorr4x4,AutocorrSum];
    
    % Window of average step (based on VT acceleration)             % we defined a window, which will help us focusing on a section of the signal to refine peaks finding
    autoCov_half = Autocorr4x4(:,1);                                % we only select the VT autocorrelation to define the window
    autoCov_half_Detr = detrend(autoCov_half);
    autoCov_half_Detr_HPFilt = HPFilter(autoCov_half_Detr,fs,Settings.Fc_HighPassFilter_AutoCorr,2);
    L = size(autoCov_half_Detr_HPFilt,1);                           % length of signal
    NFFT = 2^nextpow2(L);                                           % next power of 2 from length of signal
    Y = fft(autoCov_half_Detr_HPFilt,NFFT)/L;                       % Fourier Transform
    f = fs/2*linspace(0,1,NFFT/2+1); warning('off');                % frequencies of the signal
    Power_Spectrum = 2*abs(Y(1:NFFT/2+1));                          % single-Sided Amplitude Spectrum
    [~, Index_Dominantfrequency] = max(Power_Spectrum);             % max value of Power Spectrum corresponding to the dominant frequency and index
    Dominantfrequency = f(Index_Dominantfrequency);                 % we identify dominant frequency, to know the main period of the signal
    Window = round((1/Dominantfrequency)*fs);                       % window size will be based on the dominant cycle (step frequency)
    
    % Step & Stride Cycles
    RangeStep_Low = ceil(Window*0.5);                               % we set the lowest range (x-axes) from which we will look for the peak of step
    RangeStep_High = ceil(Window*1.5);                              % we set the highest range (x-axes) within which we will look for peak of step / this will be the lowest range to look for stride peak
    RangeStride_High = ceil(Window*2.5);                            % we set the highest range (x-axes) within which we will look for peak of sride
    
    for iChannel = 1:size(AutocorrAll_,2)-1                         % we only set our loop for 4 axes: VT, ML, AP and resultant (here, we do not look at the sum of the other 3 autocorrelations)
        AutocorrNorm = AutocorrAll_(:,iChannel)./AutocorrAll_(1,iChannel); %we normalize the signal to the value at lag = 0 (which supposes to correspond to the maximal similarity of the signal with itself)
        if iChannel == 2 %ML
            [AmplitudeAutocorr_Step_,Index_StepTimeAutocorr] = min(AutocorrNorm(RangeStep_Low:RangeStep_High,1));
            AmplitudeAutocorr_Step(:,iChannel) = abs(AmplitudeAutocorr_Step_);
            [AmplitudeAutocorr_Stride(:,iChannel),Index_StrideTimeAutocorr] = max(AutocorrNorm(RangeStep_High:RangeStride_High,1));
        else %VT, AP, Resultant
            [AmplitudeAutocorr_Step(:,iChannel),Index_StepTimeAutocorr] = max(AutocorrNorm(RangeStep_Low:RangeStep_High,1));
            [AmplitudeAutocorr_Stride(:,iChannel),Index_StrideTimeAutocorr] = max(AutocorrNorm(RangeStep_High:RangeStride_High,1));
        end
        StepTimeSamples(:,iChannel) = Index_StepTimeAutocorr + RangeStep_Low-1;
        StrideTimeSamples(:,iChannel) = Index_StrideTimeAutocorr + RangeStep_High-1;
        if iChannel == 2 %ML
            AutocorrNorm(AutocorrNorm <= 0) = AutocorrNorm(AutocorrNorm <= 0).*-1; % we transform negative part of signal into positive, so the inverted peak corresponding to the step, would give us an absolute value, to calcuate the ratio
        else %VT, AP, Resultant
            Offset = 5;                                             % offset should be constant among subjects, to be able to compare between them
            AutocorrNorm = AutocorrNorm + Offset;                   % we added an offset, to have all peaks in positive values, so the final ratio would not be biased by the initial sign of the peaks
        end
        Ad1 = AutocorrNorm(StepTimeSamples(:,iChannel),1);
        Ad2 = AutocorrNorm(StrideTimeSamples(:,iChannel),1);
        Symmetry_RatioAd1Ad2(:,iChannel) = Ad1/Ad2;                 % we calculate the ratio between the amplitude of the peaks (step/stride)
        Symmetry_DiffAd1Ad2(:,iChannel) = abs(Ad2-Ad1);             % we calculate the absolute difference between the amplitude of the peaks (step-stride)
    end
    
    %% Regularity measures
    NameChannels = {'VT','ML','AP','Resultant'};
    for iChannel = 1:size(AutocorrAll_,2)-1
        StepRegularity.(NameChannels{iChannel}) = AmplitudeAutocorr_Step(:,iChannel); % Moe-Nilssen&Helbostatt,2004
        StrideRegularity.(NameChannels{iChannel}) = AmplitudeAutocorr_Stride(:,iChannel); % Moe-Nilssen&Helbostatt,2004
        Symmetry_AutocorrelationRatio.(NameChannels{iChannel}) = Symmetry_RatioAd1Ad2(:,iChannel);
        Symmetry_AutocorrelationDiff.(NameChannels{iChannel}) = Symmetry_DiffAd1Ad2(:,iChannel);
    end
else
    for iChannel = size(NameChannels,2)
        StepRegularity.(NameChannels{iChannel}) = nan;
        StrideRegularity.(NameChannels{iChannel}) = nan;
        Symmetry_AutocorrelationRatio.(NameChannels{iChannel}) = nan;
        Symmetry_AutocorrelationDiff.(NameChannels{iChannel}) = nan;
        
    end
end
%end

% Outputs
for iChannel = [1 3]                                                       % we will focus only on VT and AP axes
    WB_SDMO.("StepRegularity_Acc_"+(NameChannels{iChannel})) = mean(extractfield(StepRegularity,(NameChannels{iChannel})));
    WB_SDMO.("StrideRegularity_Acc_"+(NameChannels{iChannel})) = mean(extractfield(StrideRegularity,(NameChannels{iChannel})));
    WB_SDMO.("Symmetry_AutocorrelationRatio_"+(NameChannels{iChannel})) = mean(extractfield(Symmetry_AutocorrelationRatio,(NameChannels{iChannel})));
    WB_SDMO.("Symmetry_AutocorrelationDiff_"+(NameChannels{iChannel})) = mean(extractfield( Symmetry_AutocorrelationDiff,(NameChannels{iChannel})));
end
end

