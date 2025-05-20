function [WB_SDMO] = extract_StepLengthVelocityVariability_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_StepLengthVelocityVariability_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of variability of step and stride length
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
%       Settings.MinRangeStepLength = 0.23;         %CONSENSUS!             % minimal range for step length
% 6) Settings for Plots created (for the development/checking of algorithm)
%       PlotSettings.StepDetectionDevelopment = 0;
%       PlotSettings.StepLengthDevelopment = 0;

%% Output
% WB_SDMO(i).StepLength_CoeffVariation = coefficient of variation of steps length for all steps within the WB (walking bout in the given imu data).
% WB_SDMO(i).StepLength_SD = standard deviation of steps length for all steps within the WB (walking bout in the given imu data).
% WB_SDMO(i).StepVelocity_CoeffVariation = coefficient of variation of  steps velocity for all steps within the WB (walking bout in the given imu data).
% WB_SDMO(i).StepVelocity_SD = standard deviation of steps velocity for all steps within the WB (walking bout in the given imu data).
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
StepLength_CoeffVariation = [];
StepLength_SD = [];
StepVelocity_CoeffVariation = [];
StepVelocity_SD = [];
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
%else                                                                % provided as an input 
%    IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
%    StepDuration = diff(IC)'./fs;                                   % [s]
%    StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);            
%end

% Spatial features calculation
SpatialFeatures = calculate_SpatialGaitFeatures_Zijlstra(AccVT,fs,IC,StepDuration,Settings,PlotSettings);
StepLength = SpatialFeatures.StepLength;
StepVelocity = SpatialFeatures.StepVelocity;
Sections= ["First", "Middle", "Last"];
% Step length & step velocity
if Settings.steps== 0
    if  StepsNumber > Settings.MinNumEvents
        Index_ExcStepLength = find(StepLength < Settings.MinRangeStepLength | StepLength > Settings.MaxRangeStepLength); % exclude anomoly step lengths
        StepLength(Index_ExcStepLength) = [];  
        for iSec= 1:length(Sections) 
            if Sections(iSec)== "First"
                First.StepLength_CoeffVariation = nanstd(StepLength(1:3))/nanmean(StepLength(1:3));
                First.StepLength_SD = nanstd(StepLength(1:3));
                %StepVelocity(Index_ExcStepLength) = [];            
                First.StepVelocity_CoeffVariation= nanstd(StepVelocity(1:3))/nanmean(StepVelocity(1:3));
                First.StepVelocity_SD = nanstd(StepVelocity(1:3));
            elseif Sections(iSec)== "Middle"
                Start=round(length(StepLength)/2)-1; End=round(length(StepLength)/2)+1;
                Middle.StepLength_CoeffVariation = nanstd(StepLength(Start:End))/nanmean(StepLength(Start:End));
                Middle.StepLength_SD = nanstd(StepLength(Start:End));
                %StepVelocity(Index_ExcStepLength) = [];  
                Start=round(length(StepVelocity)/2)-1; End=round(length(StepVelocity)/2)+1;
                Middle.StepVelocity_CoeffVariation= nanstd(StepVelocity(Start:End))/nanmean(StepVelocity(Start:End));
                Middle.StepVelocity_SD = nanstd(StepVelocity(Start:End));
            elseif Sections(iSec)== "Last"
                Last.StepLength_CoeffVariation = nanstd(StepLength(end-2:end))/nanmean(StepLength(end-2:end));
                Last.StepLength_SD = nanstd(StepLength(end-2:end));
                %StepVelocity(Index_ExcStepLength) = [];            
                Last.StepVelocity_CoeffVariation= nanstd(StepVelocity(end-2:end))/nanmean(StepVelocity(end-2:end));
                Last.StepVelocity_SD = nanstd(StepVelocity(end-2:end));
            end
        end
    else %not enough steps to define average values
        StepLength_CoeffVariation = nan;
        StepLength_SD = nan;
        StepVelocity_CoeffVariation = nan;
        StepVelocity_SD = nan;
    end
else
    if  StepsNumber > Settings.MinNumEvents
        Index_ExcStepLength = find(StepLength < Settings.MinRangeStepLength | StepLength > Settings.MaxRangeStepLength); % exclude anomoly step lengths
        StepLength(Index_ExcStepLength) = [];  
        StepLength_CoeffVariation = nanstd(StepLength)/nanmean(StepLength);
        StepLength_SD = nanstd(StepLength);
        StepVelocity(Index_ExcStepLength) = [];            
        StepVelocity_CoeffVariation= nanstd(StepVelocity)/nanmean(StepVelocity);
        StepVelocity_SD = nanstd(StepVelocity); 
    else %not enough steps to define average values
        StepLength_CoeffVariation = nan;
        StepLength_SD = nan;
        StepVelocity_CoeffVariation = nan;
        StepVelocity_SD = nan;
    end %if StepsNumber
end %Settings
%end

% Outputs
if Settings.steps== 0
    Variables= ["StepLength_CoeffVariation","StepLength_SD", "StepVelocity_CoeffVariation","StepVelocity_SD"];
    for iVar=1:length(Variables)
        for iSec = 1:length(Sections)
            VariableName= Variables(iVar);
            %VariableVector= eval(Sections(iSec).VariableName);
            %SizeVariable= length(VariableVector);
            if Sections(iSec)== "First"
                WB_SDMO.(Sections(iSec)).(VariableName) = (First.(VariableName));            
            elseif Sections(iSec)== "Middle"
                WB_SDMO.(Sections(iSec)).(VariableName)= (Middle.(VariableName));
            elseif Sections(iSec)== "Last"
                WB_SDMO.(Sections(iSec)).(VariableName)= (Last.(VariableName));
            end
     
            %WB_SDMO.(Sections(iSec)).(VariableName) = First;
            
        end
    end

else
    WB_SDMO.StepLength_CoeffVariation = nanmean(StepLength_CoeffVariation);
    WB_SDMO.StepLength_SD = nanmean(StepLength_SD);
    WB_SDMO.StepVelocity_CoeffVariation = nanmean(StepVelocity_CoeffVariation);
    WB_SDMO.StepVelocity_SD = nanmean(StepVelocity_SD);
end

