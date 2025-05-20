function [WB_SDMO] = extract_StepLengthVelocityAsymmetry_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_StepLengthVelocityAsymmetry_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of asymmetry in step length (and step velocity)
% * Ref 1) According to: McCamley John, Donati Marco. Grimpampi Eleni, Mazzà Claudia.
% "An enhanced estimate of initial contact and final contact instants of time using lower trunk inertial sensor data."
% Gait and Posture 36,(2012): 316-318"
% * Ref 2) According to: Zijlstra, Wiebren.
% "Assessment of spatio-temporal parameters during unconstrained walking."
% European journal of applied physiology 92,(2004): 39-44"
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
%       Settings.MinRangeStepLength = 0.23;         %CONSENSUS!             % minimmal range for step length
% 6) Settings for Plots created (for the development/checking of algorithm)
%       PlotSettings.StepDetectionDevelopment = 0;
%       PlotSettings.StepLengthDevelopment = 0;

%% Output
% WB_SDMO(i).StepLength_Asymmetry = asymmetry right-left length for all steps within the WB (walking bout in the given imu data).
% WB_SDMO(i).StepVelocity_Asymmetry = asymmetry right-left velocity for all steps within the WB (walking bout in the given imu data).
% Remark: Length-based features calculated from straight line walks.

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Micó Amigo, based on UNEW algorithm (BAM department, original authors: Alan Godfrey & Silvia del Din)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
       
%% History
% 2019/6th/September functionized - Encarna Micó Amigo


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
StepLength_Asymmetry = [];
StepVelocity_Asymmetry = [];
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB

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

% Spatial features calculation
SpatialFeatures = calculate_SpatialGaitFeatures_Zijlstra(AccVT,fs,IC,StepDuration,Settings,PlotSettings);
StepLength = SpatialFeatures.StepLength;
StepVelocity = SpatialFeatures.StepVelocity;

% Step length & step velocity
if  StepsNumber > Settings.MinNumEvents
    Index_ExcStepLength = find(StepLength < Settings.MinRangeStepLength | StepLength > Settings.MaxRangeStepLength); % exclude anomoly step lengths
    % Step length
    StepLength(Index_ExcStepLength) = [];            
    StepLength_R = StepLength(1:2:end);
    StepLength_L = StepLength(2:2:end);
    if (size(StepLength_R,1) >= 1) && (size(StepLength_L,1) >= 1)
        StepLength_Asymmetry = abs(mean(StepLength_R) - mean(StepLength_L));
    end
    % Step velocity
    StepVelocity(Index_ExcStepLength) = [];
    StepVelocity_R = StepVelocity(1:2:end);
    StepVelocity_L = StepVelocity(2:2:end);
    if (size(StepVelocity_R,1) >= 1) && (size(StepVelocity_L,1) >= 1)
        StepVelocity_Asymmetry = abs(mean(StepVelocity_R) - mean(StepVelocity_L));
    end
else %not enough steps to define average values
    StepLength_Asymmetry = nan;
    StepVelocity_Asymmetry = nan;
end
%end

% Outputs
WB_SDMO.StepLength_Asymmetry = nanmean(StepLength_Asymmetry);
WB_SDMO.StepVelocity_Asymmetry = nanmean(StepVelocity_Asymmetry);
end

