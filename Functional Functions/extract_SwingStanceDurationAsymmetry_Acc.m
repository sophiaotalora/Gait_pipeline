function [WB_SDMO] = extract_SwingStanceDurationAsymmetry_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_SwingStanceDurationAsymmetry_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of swing and stance phases asymmetry
% * Ref 1) According to: McCamley John, Donati Marco. Grimpampi Eleni, Mazz� Claudia.
% "An enhanced estimate of initial contact and final contact instants of time using lower trunk inertial sensor data."
% Gait and Posture 36,(2012): 316-318"
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
%       Settings.GravityConstant = 9.81; % [m/s2]                                           
%       Settings.Fc_LowPassFilter = 20*fs/100;                              % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.Fc_HighPassFilter = 0.1*fs/100;                            % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.Fc_HighPassFilter_AutoCorr = 0.8*fs/100;                   % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
%       Settings.StepDetectMcCameley_EMA = 1;                               % use McCamely algorithm with Encarna modifications (1) or with BAM original (0)
%       Settings for step length calculation
%       Settings.L5_height = 1.05; % [m]                                    % height of the sensor: from the floor to L5
%       Settings.CorrectionFactor = 1.25;                                   % correction factor defined by Zijlstra, based no the iverted pendulum model
%       Settings.MaxRangeStepLength = 0.85;         %CONSENSUS!             % maximal range for step length
%       Settings.MinRangeStepLength = 0.23;         %CONSENSUS!             % minimal range for step length
%       Settings.MaxRangeStepDuration = 1.25;       %CONSENSUS!             % maximal range for step duration
%       Settings.MinRangeStepDuration = 0.25;       %CONSENSUS!             % minimal range for step duration
%       Settings.MaxRangeSwingDuration = 0.85;      %CONSENSUS!             % maximal range for swing duration
%       Settings.MinRangeSwingDuration = 0.23;      %CONSENSUS!             % minimal range for swing duration 
% 6) Settings for Plots created (for the development/checking of algorithm)
%       PlotSettings.StepDetectionDevelopment = 0;
%       PlotSettings.StepLengthDevelopment = 0;

%% Output
% WB_SDMO(i).SwingDuration_Asymmetry = asymmetry right-left for all swing durations within the WB (walking bout in the given imu data).
% WB_SDMO(i).StanceDuration_Asymmetry = asymmetry right-left for all stance durations within the WB (walking bout in the given imu data).

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Mic� Amigo, based on UNEW algorithm (BAM department, original authors: Alan Godfrey & Silvia del Din)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
       
%% History
% 2019/6th/September functionized - Encarna Mic� Amigo


 %% SETTINGS 

 if nargin < 5                                                              % if the settings are not predefined:
     % Settings for Spatiotemporal extraction
     Settings.MinNumEvents = 3;                  %CONSENSUS!                % minimal number of steps to consider the step-to-stpe gait features
     Settings.PercOfDominantFreq = 0.5;                                     % percentage of the dominant frequency considered as a threshold (minimal distance between peaks)
     Settings.GravityConstant = 9.81; % [m/s2]
     Settings.Fc_LowPassFilter = 20*fs/100;                                 % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
     Settings.Fc_HighPassFilter = 0.1*fs/100;                               % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
     Settings.Fc_HighPassFilter_AutoCorr = 0.8*fs/100;                      % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
     Settings.StepDetectMcCameley_EMA = 1;                                  % use McCamely algorithm with Encarna modifications (1) or with BAM original (0)
     % Settings for step length calculation
     Settings.L5_height = 1.05; % [m]                                       % height of the sensor: from the floor to L5
     Settings.CorrectionFactor = 1.25;                                      % correction factor defined by Zijlstra, based no the iverted pendulum model
     Settings.MaxRangeStepLength = 0.85;         %CONSENSUS!                % maximal range for step length
     Settings.MinRangeStepLength = 0.23;         %CONSENSUS!                % minimal range for step length
     Settings.MaxRangeStepDuration = 1.25;       %CONSENSUS!                % maximal range for step duration
     Settings.MinRangeStepDuration = 0.25;       %CONSENSUS!                % minimal range for step duration
     Settings.MaxRangeSwingDuration = 0.85;      %CONSENSUS!                % maximal range for swing duration
     Settings.MinRangeSwingDuration = 0.23;      %CONSENSUS!                % minimal range for swing duration
     % PlotSettings
     PlotSettings.StepDetectionDevelopment = 0;
     PlotSettings.StepLengthDevelopment = 0;
 end

%% LOOP OVER WALKING BOUTS

%WB_Number = size(WB_WBD,2);                                               %number of walking bouts for loop
%for WB_index = 1:WB_Number                                                %loop over each WB
    
% Info about straight line walks within WB
%StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);              %number/amount of straight walk episodes
straightWalk_start = cell2mat({WB_WBD.start});                             %start of straight WB
straightWalk_end = cell2mat({WB_WBD.end});                                 %end of straight WB
    
% Initialization of outcomes
SwingDuration_Asymmetry = [];
StanceDuration_Asymmetry = [];
%     for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB

% Definition of signal: start-end, the sign of the signal does not affect the result
AccVT = imu.acc_V(straightWalk_start:straightWalk_end);
        
% Initial temporal features calculation
%if ~isfield(WB_WBD,'TempEvents')                                    % not provided as an input (FAU & EPFL work)
if Settings.StepDetectMcCameley_EMA
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
else
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
end
IC = TemporalFeatures.InitialContactSample;
StepDuration = TemporalFeatures.StepDuration;
StepsNumber = size(IC,2); 
%         else                                                                % provided as an input 
%             IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
%             StepDuration = diff(IC)'./fs;                                   % [s]
%             StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);            
%         end

% Temporal features
SwingDuration = TemporalFeatures.SwingDuration;
StanceDuration = TemporalFeatures.StanceDuration;

% Swing phase duration & stance phase duration
if  StepsNumber > Settings.MinNumEvents
    % Swing phase durations
    Index_ExcSwingDuration = find(SwingDuration < Settings.MinRangeSwingDuration | SwingDuration > Settings.MaxRangeSwingDuration); % exclude anomoly swing durations
    SwingDuration(Index_ExcSwingDuration) = [];
    SwingDuration_R = SwingDuration(1:2:end);
    SwingDuration_L = SwingDuration(2:2:end);
    if (size(SwingDuration_R,1) >= 1) && (size(SwingDuration_L,1) >= 1)
        SwingDuration_Asymmetry = abs(mean(SwingDuration_R) - mean(SwingDuration_L));
    end 
    % Stance phase durations
    StanceDuration(Index_ExcSwingDuration) = [];
    StanceDuration_R = StanceDuration(1:2:end);
    StanceDuration_L = StanceDuration(2:2:end);
    if (size(StanceDuration_R,1) >= 1) && (size(StanceDuration_L,1) >= 1)
        StanceDuration_Asymmetry = abs(mean(StanceDuration_R) - mean(StanceDuration_L));
    end
else %not enough steps to define mean values
    SwingDuration_Asymmetry = nan;
    StanceDuration_Asymmetry = nan;
end
%end
% Outputs
WB_SDMO.SwingDuration_Asymmetry = nanmean(SwingDuration_Asymmetry);
WB_SDMO.StanceDuration_Asymmetry = nanmean(StanceDuration_Asymmetry);
end

