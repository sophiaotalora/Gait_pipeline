function [WB_SDMO] = extract_Regularity_GaitSymmetryIndex_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_Regularity_GaitSymmetryIndex_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Regularity measures based on Autocorrelation signal - ACCELEROMETRY (VT and AP)
% * Ref 1) "Gait Symmetry Assessment with a Low Back 3D Accelerometer in Post-Stroke Patients"
% Authors: Zhang W1, Smuck M2,3, Legault C4, Ith MA5,6, Muaremi A7, Aminian K8.
% 2019, Journal "Sensors". https://www.ncbi.nlm.nih.gov/pubmed/30282947
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
% WB_SDMO(i).GaitSymmetryIndex_Autocorr = Gait Symmetry Index (new approach based on the combination of accelerometry)

%% Remarks
% *** Comments located on the right side
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

%WB_Number = size(WB_WBD,2);                                               % number of walking bouts for loop
%for WB_index = 1:WB_Number                                                % loop over each WB
    
% Info about straight line walks within WB
%StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);              % number/amount of straight walk episodes
straightWalk_start = cell2mat({WB_WBD.start});                             % start of straight WB
straightWalk_end = cell2mat({WB_WBD.end});                                 %end of straight WB
    
% Initialization of outcomes
GaitSymmetryIndex_Autocorr = [];
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
    
% Definition of signal: start-end
Acc = zeros(size(imu.acc_V(straightWalk_start:straightWalk_end),1),4);
Acc(:,1) = imu.acc_V(straightWalk_start:straightWalk_end);
Acc(:,2) = imu.acc_ML(straightWalk_start:straightWalk_end);
Acc(:,3) = imu.acc_AP(straightWalk_start:straightWalk_end);
if mean(Acc(:,1)) < 0                                               % signal negative = VT acceleration input is downwards
    Acc(:,1) = -Acc(:,1);                                           % we invert the signal (to get an upwards direction). If it would be directly provided as upwards, there would not be need to invert it
    Acc(:,3) = -Acc(:,3);                                           % we invert the signal (to get an forwards direction). If it would be directly provided as forwards, there would not be need to invert it
end %although there is no need no invert it

% Detrend & filter signals *** Fc_LowPassFilter is 1/2
Acc_Detr = detrend(Acc);                                            % remove the trend/offset of gravity
Acc_Detr_LPFilt = WintFilt_low(Acc_Detr-mean(Acc_Detr),Settings.Fc_LowPassFilter/2, fs); %low-pass filter the signal
Acc_Detr_LPFilt(:,4) = sqrt(sum(Acc_Detr_LPFilt(:,1:3)'.^2)');      % resultant or combined signal

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

% Calculate Regularity based on Autocorrelation signal - Acc
NameChannels = {'VT','ML','AP','Combined'};       

% Calculate Power Spectral Density - Acc
if StepsNumber > Settings.MinNumEvents
    % Stride time and regularity from auto correlation
    [Autocorr4x4,Lags] = xcov(Acc_Detr_LPFilt,'biased');            % auto-covariance of Moe-Nilssen. "Biased": attenuated signal 
    Autocorr4x4 = Autocorr4x4(Lags >= 0,[1 6 11 16]);               % we get positive lags signal (Lags>=0), we select channels 1, 6, 11 and 16 because the function has as an output the correlation of each of the signals with each of the signals, thus, to get the "auto"correlation, we select particular axes of the final matrix
    Autocorr4x4_Norm = Autocorr4x4./Autocorr4x4(1,:);
    % Set the size of the signals, take up to the first 400 samples
    if size(Autocorr4x4_Norm,1) < 400
        Range = 1:size(Autocorr4x4_Norm,1);
    else
        Range = 1:400; 
    end
    
    % All negative part of the signal should be 0
    Autocorr4x4_1stSection = Autocorr4x4_Norm(Range,:);
    Autocorr4x4_1stSection(Autocorr4x4_1stSection < 0) = 0;
    
    % Get signals based on axes
    AutocorrSum = sum(Autocorr4x4_1stSection(:,1:3),2); %Cstride    % sum the first axes-signals of autocorrelation
    AutocorrSquared = sqrt(sum(Autocorr4x4_1stSection(:,1:3),2)); %Cstep   % square of the sum
    
    % Identify the samples at which there is the maximas over signals
    [TStrideValue, TStrideSample] = max(AutocorrSum(10:end,1));     % we search 10 samples later, to avoid getting the first peak at lag = 0
    TStride = TStrideSample +10-1;
    GaitSymmetryIndex_Autocorr = AutocorrSquared(ceil(0.5*TStride),1)./sqrt(3);
else
    GaitSymmetryIndex_Autocorr = nan;
end
%end

% Outputs
WB_SDMO = mean(GaitSymmetryIndex_Autocorr);
end

