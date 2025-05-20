function [WB_SDMO] = extract_StepLengthVelocityMean_Acc(imu,fs,WB,WB_SDMO, Settings,PlotSettings)
%[WB_SDMO] = extract_StepLengthVelocityMean_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of length of step and stride cycles
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
%       Settings.L5_height = 1.05; % [m]                                    % height of the sensor: from the floor to L5
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
% 6) Settings for Plots created (for the development/checking of algorithm)
%       PlotSettings.StepDetectionDevelopment = 0;
%       PlotSettings.StepLengthDevelopment = 0;

%% Output
% WB_SDMO(i).StepLength_Mean = mean step length for all steps within the WB (walking bout in the given imu data).
% WB_SDMO(i).StepVelocity_Mean = mean step velocity for all steps within the WB (walking bout in the given imu data).
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

%WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
%for WB_index = 1:WB_Number                                                  % loop over each WB
    
% Info about straight line walks within WB
%StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);            % number/amount of straight walk episodes
straightWalk_start = cell2mat({WB.start});   % start of straight WB
straightWalk_end = cell2mat({WB.end});       % end of straight WB

% Initialization of outcomes
StepLength_Mean = [];
StepVelocity_Mean = [];
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB

% Definition of signal: start-end, the sign of the signal does not affect the result
AccVT = imu.acc_V(straightWalk_start:straightWalk_end);
        
% Initial temporal features calculation
%if ~isfield(WB,'TempEvents')                                    % not provided as an input (FAU & EPFL work)
if Settings.StepDetectMcCameley_EMA == 1
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
elseif Settings.StepDetectMcCameley_EMA == 0
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
    IC = TemporalFeatures.InitialContactSample;
    StepDuration = TemporalFeatures.StepDuration;
    StepsNumber = size(IC,2); 
elseif Settings.StepDetectMcCameley_EMA == 3
    IC= round(straightWalk_start/fs);
    StepDuration = (straightWalk_end/fs) - (straightWalk_start/fs);
    StepsNumber = 1;
end
%IC = TemporalFeatures.InitialContactSample;
%StepDuration = TemporalFeatures.StepDuration;
%StepsNumber = size(IC,2); 
%else                                                                % provided as an input 
%    IC = WB(WB_index).straightWalk.TempEvents;                  % [sample number]
%    StepDuration = diff(IC)'./fs;                                   % [s]
%    StepsNumber = WB(WB_index).straightWalk.stepsNumber;        % size(IC,2);            
%end

% Spatial features calculation
SpatialFeatures = calculate_SpatialGaitFeatures_Zijlstra(AccVT,fs,IC,StepDuration,Settings,PlotSettings);
StepLength = SpatialFeatures.StepLength;
StepVelocity = SpatialFeatures.StepVelocity;

% Step length & step velocity
if  StepsNumber > Settings.MinNumEvents 
    Index_ExcStepLength = find(StepLength < Settings.MinRangeStepLength | StepLength > Settings.MaxRangeStepLength); % exclude anomoly step lengths
    StepLength(Index_ExcStepLength) = [];
    StepLength_Mean = nanmean(StepLength);
    StepVelocity(Index_ExcStepLength) = [];
    StepVelocity_Mean = nanmean(StepVelocity);
else %not enough steps to define mean values
    StepLength_Mean = nan;
    StepVelocity_Mean = nan;
end
%end
% Outputs
if Settings.steps== 0
    Variables= ["StepLength","StepVelocity"];
    Sections= ["First", "Middle", "Last"];
    for iVar=1:length(Variables)
        VariableName= Variables(iVar);
        VariableVector= eval(VariableName);
        SizeVariable= length(VariableVector);
        if contains(VariableName, "Stride") || contains(VariableName, "stride")
            First = VariableVector(1:2);
            Middle = VariableVector(round(SizeVariable/2):round(SizeVariable/2)+1);
            Last = VariableVector(end-1:end);      
        else
            First = VariableVector(1:3);
            Middle = VariableVector(round(SizeVariable/2)-1:round(SizeVariable/2)+1);
            Last = VariableVector(end-2:end);
        end

        WB_SDMO.(Sections(1)).(VariableName) = First;
        WB_SDMO.(Sections(2)).(VariableName) = Middle;
        WB_SDMO.(Sections(3)).(VariableName) = Last;
    end

else
    WB_SDMO.StepLength_Mean = mean(StepLength_Mean);
    WB_SDMO.StepVelocity_Mean = mean(StepVelocity_Mean);
end
% Outputs
%WB_SDMO.StepLength_Mean = mean(StepLength_Mean);
%WB_SDMO.StepVelocity_Mean = mean(StepVelocity_Mean);
%end
end
