function [WB_SDMO] = extract_LyapunovExponent_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_LyapunovExponent_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Lyapunov exponent - ANGULAR VELOCITY / GYROSCOPE
% * Ref 1) Authors: Arnaud Dupeyron,Sietse M. Rispens, Christophe Demattei and Jaap H. van Dieën
% "Precision of estimates of local stability of repetitive trunk movements"
% Eur Spine J. 2013 Dec; 22(12): 2678–2685.
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
%       Settings.PercOfDominantFreq = 0.5;                                  % percentage of the dominant frequency considered as a threshold (minimal distance between peaks)
%       Settings.L5_height = 1.05; % [m]                                    % height of the sensor: from the floor to L5
%       Settings.GravityConstant = 9.81; % [m/s2]
%       Settings.Fc_LowPassFilter = 20*fs/100;                              % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.Fc_HighPassFilter = 0.1*fs/100;                            % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.Fc_HighPassFilter_AutoCorr = 0.8*fs/100;                   % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.StepDetectMcCameley_EMA = 1;                               % use McCamely algorithm with Encarna modifications (1) or with BAM original (0)
% 6) Settings for phase/complex measures
%       Settings.MinNumStrides_Complex = round(Settings.MinNumEvents/2);    % minimal number of strides recommended
% 7) Settings for Plots created (for the development/checking of algorithm)
%       PlotSettings.PhaseFeatures = 0;

%% Output
% WB_SDMO(i).LyapunovRC_Gyr = Local dynamic stability exponent (Rosestein method) from angular velocity data, for each signal-axis 
% WB_SDMO(i).LyapunovW_Gyr = Local dynamic stability exponent (Wolf method) from angular velocity data, for each signal-axis 

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
    % Settings for phase/complex measures
    Settings.MinNumStrides_Complex = floor(Settings.MinNumEvents/2);        % minimal number of strides in the signal
    % PlotSettings
    PlotSettings.PhaseFeatures = 0;                                         % activate plots to help understangind the development of the phase analysis algorithms
end

%% LOOP OVER WALKING BOUTS

%WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
%for WB_index = 1:WB_Number                                                  % loop over each WB
    
% Info about straight line walks within WB
%StraightWalk_Number = size(WB_WBD.straightWalk,2);            % number/amount of straight walk episodes
straightWalk_start = cell2mat({WB_WBD.start});   % start of straight WB
straightWalk_end = cell2mat({WB_WBD.end});       % end of straight WB
NameChannels = {'Yaw','Pitch','Roll','Combined'};

% Initialization
clear LyapunovRC
clear LyapunovW
for iChannel = 1:size(NameChannels,2)
    LyapunovRC.(NameChannels{iChannel}) = [];
    LyapunovW.(NameChannels{iChannel}) = [];
end
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
    
% Definition of signal: start-end
Gyr = zeros(size(imu.gyro_yaw(straightWalk_start:straightWalk_end),1),4);
Gyr(:,1) = imu.gyro_yaw(straightWalk_start:straightWalk_end);
Gyr(:,2) = imu.gyro_pitch(straightWalk_start:straightWalk_end);
Gyr(:,3) = imu.gyro_roll(straightWalk_start:straightWalk_end);

% Detrend & filter signals
Gyr_Detr = detrend(Gyr);                                            % remove the trend
Gyr_Detr_LPFilt = WintFilt_low(Gyr_Detr-mean(Gyr_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal
Gyr_Detr_LPFilt(:,4) = sqrt(sum(Gyr_Detr_LPFilt(:,1:3)'.^2)');      % resultant or combined signal

% Initial temporal features calculation
% if ~isfield(WB_WBD,'TempEvents')                                     % not provided as an input (FAU & EPFL work)
AccVT = imu.acc_V(straightWalk_start:straightWalk_end);
if Settings.StepDetectMcCameley_EMA
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
else
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
end
IC = TemporalFeatures.InitialContactSample;
StepDuration = TemporalFeatures.StepDuration;
StepsNumber = size(IC,2);
% else                                                                % provided as an input
%     IC = WB_WBD.straightWalk.TempEvents;                  % [sample number]
%     StepDuration = diff(IC)'./fs;                                   % [s]
%     StepsNumber = WB_WBD.straightWalk.stepsNumber;        % size(IC,2);
% end

% Calculation of phase features calculation - Acc
if  StepsNumber > Settings.MinNumStrides_Complex
    % These settings can be pre-calculated with appropriate (also available) functions. However, for the time being, we pre-set them as follows:
    Ly_J = round(12*fs/100); % samples                              % Embedding delay (used in Lyapunov estimations)
    Ly_m = 5; % dimensions                                          % Embedding dimension (used in Lyapunov estimations)
    Ly_FitWinLen = round(60*fs/100); %Window                        % Fitting window length (used in Lyapunov estimations Rosenstein's method)
    for iChannel = 1:size(Gyr_Detr_LPFilt,2)                
        % Rosestein
        [LyapunovRC.(NameChannels{iChannel}),~] = CalcMaxLyapConvGait(Gyr_Detr_LPFilt(:,iChannel),fs,struct('J',Ly_J,'m',Ly_m,'FitWinLen',Ly_FitWinLen));
        % Wolf
        [LyapunovW.(NameChannels{iChannel}),~] = CalcMaxLyapWolfFixedEvolv(Gyr_Detr_LPFilt(:,iChannel),fs,struct('J',Ly_J,'m',Ly_m));
    end
else %not enough steps to define average values
    for iChannel = 1:size(Gyr_Detr_LPFilt,2)
        LyapunovRC.(NameChannels{iChannel}) = nan;
        LyapunovW.(NameChannels{iChannel}) = nan;
    end
end
%end

% Outputs
for iChannel = 1:size(Gyr_Detr_LPFilt,2)
    WB_SDMO.LyapunovRC_Gyr.(NameChannels{iChannel}) = nanmean(extractfield(LyapunovRC,(NameChannels{iChannel})));
    WB_SDMO.LyapunovW_Gyr.(NameChannels{iChannel}) = nanmean(extractfield(LyapunovW,(NameChannels{iChannel})));
end
%end

