function [WB_SDMO] = extract_contentSpectralPowerDensity_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_contentSpectralPowerDensity_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Spectral Components - GYROSCOPE / ANGULAR VELOCITY
% * Ref 1) Authors: Weiss A1, Sharifi S, Plotnik M, van Vugt JP, Giladi N, Hausdorff JM.
% "Toward automated, at-home assessment of mobility among patients with Parkinson disease, using a body-worn accelerometer."
% Neurorehabil Neural Repair. 2011 Nov-Dec;25(9):810-8
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
%       PlotSettings.SpectralFatures = 1;                                   % activate plots to help understangind the calculation of power spectral density of signals

%% Output
% WB_SDMO(i).SpectralFeatures_Gyro = spectral components from angular velocity signals (for each signal-axis)

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Micó Amigo, based on Vrije Universiteit Amsterdam algorithm (original author: Sietse Rispens)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com

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
    PlotSettings.SpectralFatures = 1;                                       % activate plots to help understangind the calculation of power spectral density of signals
end

%% LOOP OVER WALKING BOUTS

%WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
%for WB_index = 1:WB_Number                                                  % loop over each WB
    
% Info about straight line walks within WB
%    StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);            % number/amount of straight walk episodes
straightWalk_start = cell2mat({WB_WBD.start});   % start of straight WB
straightWalk_end = cell2mat({WB_WBD.end});       % end of straight WB
NameChannels = {'Yaw','Pitch','Roll'};

% Initialization
VariablesPhaseNames = {'DominantFreq','Amplitude','AmplitudeNorm','Width','Slope','Range','MeanPower','MedianPower','MeanFreq','MedianFreq','WidthNorm','SlopeNorm','RangeNorm','IntegratedPower','HarmonicRatio','IndexHarmonicity'};
clear SpectralFeatures
for iChannel = 1:size(NameChannels,2)
    for iFeature = 1:size(VariablesPhaseNames,1)
        SpectralFeatures.(NameChannels{iChannel}).(VariablesPhaseNames{iFeature}) = [];
    end
end
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
    
% Definition of signal: start-end
% Gyr
Gyr = zeros(size(imu.gyro_yaw(straightWalk_start:straightWalk_end),1),3);
Gyr(:,1) = imu.gyro_yaw(straightWalk_start:straightWalk_end);
Gyr(:,2) = imu.gyro_pitch(straightWalk_start:straightWalk_end);
Gyr(:,3) = imu.gyro_roll(straightWalk_start:straightWalk_end);
% Acc
Acc = zeros(size(imu.acc_V(straightWalk_start:straightWalk_end),1),3);
Acc(:,1) = imu.acc_V(straightWalk_start:straightWalk_end);
Acc(:,2) = imu.acc_ML(straightWalk_start:straightWalk_end);
Acc(:,3) = imu.acc_AP(straightWalk_start:straightWalk_end);
if mean(Acc(:,1)) < 0                                               % signal negative = VT acceleration input is downwards
    Acc(:,1) = -Acc(:,1);                                           % we invert the signal (to get an upwards direction). If it would be directly provided as upwards, there would not be need to invert it
    Acc(:,3) = -Acc(:,3);                                           % we invert the signal (to get an forwards direction). If it would be directly provided as forwards, there would not be need to invert it
end %although there is no need no invert it

% Detrend & filter signals
Gyr_Detr = detrend(Gyr);                                            % remove the trend
Gyr_Detr_LPFilt = WintFilt_low(Gyr_Detr-mean(Gyr_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal
Acc_Detr = detrend(Acc);                                            % remove the trend/offset of gravity
Acc_Detr_LPFilt = WintFilt_low(Acc_Detr-mean(Acc_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal

% Initial temporal features calculation
%if ~isfield(WB_WBD,'TempEvents')                                    % not provided as an input (FAU & EPFL work)
AccVT = Acc(:,1);
if Settings.StepDetectMcCameley_EMA
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
else
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
end
IC = TemporalFeatures.InitialContactSample;
StepsNumber = size(IC,2);
%else                                                                % provided as an input
%    IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
%    StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);
%end

% Calculate Power Spectral Density - Gyr
if StepsNumber > Settings.MinNumEvents
    %Define minimal straight walking bout duration, since we need the same junk for each signal of comparison
    %WB = [];
    %for iWB = 1:size(WB_WBD,2)
    %    StraightWB = size(WB_WBD(iWB).straightWalk,2);
    %    for iSWB = 1:StraightWB
    %        SamplesTotalSig = WB_WBD(iWB).straightWalk(iSWB).end - WB_WBD(iWB).straightWalk(iSWB).start;
    %        WB = [WB; SamplesTotalSig];
    %    end
    %end
    %HanningWindow_Size = min(WB);                                   % shortest walking bout, for comparison, all signals need to have the same length
    HanningWindow_Size= Settings.HanningSize*fs;
    Nfft = HanningWindow_Size*10;
    Signal_iD = 2;                                                  % 1 = accelerometry, 2 = gyroscope
    % get Stride Frequency (based on Accelerometry)
    for iAcc = 1:size(Acc_Detr_LPFilt,2)
        HanninWindow = hanning(floor(HanningWindow_Size));
        [PowerSpectralDensity(:,iAcc),Frequencies] = pwelch(Acc_Detr_LPFilt(1:HanningWindow_Size,iAcc),HanninWindow,50,Nfft,fs); %normalized
    end
    dF = Frequencies(2)-Frequencies(1); %resolution: frequency/bit
    StrideDurationSeconds = mean(TemporalFeatures.StrideDuration);
    [StrideFrequency, QualityIndicator, ~] = strideFrequencyFrom3dAcc(PowerSpectralDensity(:,1:3),Frequencies);
    if QualityIndicator < 1 && isfinite(StrideDurationSeconds) && StrideDurationSeconds ~= 0
        StrideFrequency = 1./StrideDurationSeconds;
    end
    % get Frequency Analysis - Gyroscope
    SpectralFeatures_Gyr = getFrequencyAnalysis...
        (Gyr_Detr_LPFilt,Gyr,fs,Nfft,Settings.N_Harm,Settings,StrideFrequency,Signal_iD,HanningWindow_Size,PlotSettings);
    FeaturesName = fieldnames(SpectralFeatures_Gyr);
    FeaturesNumMax = size(FeaturesName,1);
    for iFeature = 1:FeaturesNumMax
        for iChannel = 1:size(SpectralFeatures_Gyr.(FeaturesName{iFeature}),2)
            SpectralFeatures.(NameChannels{iChannel}).(VariablesPhaseNames{iFeature}) = ...
                SpectralFeatures_Gyr.(FeaturesName{iFeature})(iChannel);
        end
    end
else
    for iChannel = 1:size(Gyr_Detr_LPFilt,2)
        VariablesPhaseNames = fieldnames(SpectralFeatures.(NameChannels{iChannel}));
        for iFeature = 1:size(VariablesPhaseNames,1)
            SpectralFeatures.(NameChannels{iChannel}).(VariablesPhaseNames{iFeature}) = nan;
        end
    end
end
%end
    
% Outputs
VariablesPhaseNames = fieldnames(SpectralFeatures.(NameChannels{iChannel})); %it will take the variable names (list) from the last straight walk bout and the last channel
for iChannel = 1:size(Gyr_Detr_LPFilt,2)
    for iFeature = 1:size(VariablesPhaseNames,1)
        %DataStraightWalk = [];
        %for iStraightWalks = 1:StraightWalk_Number %size(SpectralFeatures,2) % we need to accumulate all variables for each straight walk within the WB, and then calculate the mean value
        DataStraightWalk_ = SpectralFeatures.(NameChannels{iChannel}).(VariablesPhaseNames{iFeature}');
        %    DataStraightWalk = [DataStraightWalk; DataStraightWalk_];
        %end
        WB_SDMO.SpectralFeatures_Gyr.(NameChannels{iChannel}).(VariablesPhaseNames{iFeature}) = mean(DataStraightWalk_);
        %clear DataStraightWalk_
    end
end
end


