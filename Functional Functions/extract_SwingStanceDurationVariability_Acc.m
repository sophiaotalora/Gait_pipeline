function [WB_SDMO] = extract_SwingStanceDurationVariability_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_SwingStanceDurationVariability_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of variability of swing and stance cycles duration
% * Ref 1) According to: McCamley John, Donati Marco. Grimpampi Eleni, Mazzà Claudia.
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
% WB_SDMO(i).SwingDuration_CoeffVariation = coefficient of variation of swing durations for all strides within the WB (walking bout in the given imu data).
% WB_SDMO(i).SwingDuration_SD = standard deviation of swing durations for all strides within the WB (walking bout in the given imu data).
% WB_SDMO(i).StanceDuration_CoeffVariation = coefficient of variation of stance durations for all strides within the WB (walking bout in the given imu data).
% WB_SDMO(i).StanceDuration_SD = standard deviation of stance durations for all strides within the WB (walking bout in the given imu data).

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
SwingDuration_CoeffVariation = [];
SwingDuration_SD = [];
SwingDurationPercentange_CoeffVariation = [];
SwingDurationPercentange_SD = [];
StanceDuration_CoeffVariation = [];
StanceDuration_SD = [];
StanceDurationPercentange_CoeffVariation = [];
StanceDurationPercentange_SD = [];
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

% Temporal features
SwingDuration = TemporalFeatures.SwingDuration;
StanceDuration = TemporalFeatures.StanceDuration;
SwingDurationPercentage = (SwingDuration./TemporalFeatures.StrideDuration)*100;
StanceDurationPercentage = (StanceDuration/TemporalFeatures.StrideDuration)*100;
Sections= ["First", "Middle", "Last"];

% Swing phase duration & stance phase duration
if Settings.steps== 0
    if  StepsNumber > Settings.MinNumEvents
        Index_ExcSwingDuration = find(SwingDuration < Settings.MinRangeSwingDuration | SwingDuration > Settings.MaxRangeSwingDuration); % exclude anomoly swing durations
        % Swing
        SwingDuration(Index_ExcSwingDuration) = [];
        SwingDurationPercentage(Index_ExcSwingDuration) = [];
        StanceDuration(Index_ExcSwingDuration) = [];
        StanceDurationPercentage(Index_ExcSwingDuration) = [];
        for iSec= 1:length(Sections) 
                if Sections(iSec)== "First"
                    SwingDuration_CoeffVariation = nanstd(SwingDuration(1:3))/nanmean(SwingDuration(1:3));
                    SwingDurationPercentage_CoeffVariation = nanstd(SwingDurationPercentage(1:3))/nanmean(SwingDurationPercentage(1:3));
                    SwingDuration_SD = nanstd(SwingDuration(1:3));
                    SwingDurationPercentage_SD = nanstd(SwingDurationPercentage(1:3));
                    % Stance
                    StanceDuration_CoeffVariation = nanstd(StanceDuration(1:3))/nanmean(StanceDuration(1:3));
                    StanceDurationPercentage_CoeffVariation = nanstd(StanceDurationPercentage(1:3))/nanmean(StanceDurationPercentage(1:3));
                    StanceDuration_SD = nanstd(StanceDuration(1:3));
                    StanceDurationPercentage_SD = nanstd(StanceDurationPercentage(1:3));
                elseif Sections(iSec)== "Middle" 
                    SwingDuration_CoeffVariation = nanstd(SwingDuration(1:3))/nanmean(SwingDuration(1:3));
                    SwingDurationPercentage_CoeffVariation = nanstd(SwingDurationPercentage(1:3))/nanmean(SwingDurationPercentage(1:3));
                    SwingDuration_SD = nanstd(SwingDuration(1:3));
                    SwingDurationPercentage_SD = nanstd(SwingDurationPercentage(1:3));
                    % Stance
                    StanceDuration_CoeffVariation = nanstd(StanceDuration(1:3))/nanmean(StanceDuration(1:3));
                    StanceDurationPercentage_CoeffVariation = nanstd(StanceDurationPercentage(1:3))/nanmean(StanceDurationPercentage(1:3));
                    StanceDuration_SD = nanstd(StanceDuration(1:3));
                    StanceDurationPercentage_SD = nanstd(StanceDurationPercentage(1:3));

                elseif Sections(iSec)== "Last"
                    SwingDuration_CoeffVariation = nanstd(SwingDuration(1:3))/nanmean(SwingDuration(1:3));
                    SwingDurationPercentage_CoeffVariation = nanstd(SwingDurationPercentage(1:3))/nanmean(SwingDurationPercentage(1:3));
                    SwingDuration_SD = nanstd(SwingDuration(1:3));
                    SwingDurationPercentage_SD = nanstd(SwingDurationPercentage(1:3));
                    % Stance
                    StanceDuration_CoeffVariation = nanstd(StanceDuration(1:3))/nanmean(StanceDuration(1:3));
                    StanceDurationPercentage_CoeffVariation = nanstd(StanceDurationPercentage(1:3))/nanmean(StanceDurationPercentage(1:3));
                    StanceDuration_SD = nanstd(StanceDuration(1:3));
                    StanceDurationPercentage_SD = nanstd(StanceDurationPercentage(1:3));

                end
        end
    else %not enough steps to define average values
        % Swing
        SwingDuration_CoeffVariation = nan;
        SwingDurationPercentage_CoeffVariation = nan;
        SwingDuration_SD = nan;
        SwingDurationPercentage_SD = nan;
        % Stance
        StanceDuration_CoeffVariation = nan;
        StanceDurationPercentage_CoeffVariation = nan;
        StanceDuration_SD = nan;
        StanceDurationPercentage_SD = nan;   
        
    end %if StepsNumber
else
    if  StepsNumber > Settings.MinNumEvents
        Index_ExcSwingDuration = find(SwingDuration < Settings.MinRangeSwingDuration | SwingDuration > Settings.MaxRangeSwingDuration); % exclude anomoly swing durations
        % Swing
        SwingDuration(Index_ExcSwingDuration) = [];
        SwingDurationPercentage(Index_ExcSwingDuration) = [];
        SwingDuration_CoeffVariation = nanstd(SwingDuration)/nanmean(SwingDuration);
        SwingDurationPercentage_CoeffVariation = nanstd(SwingDurationPercentage)/nanmean(SwingDurationPercentage);
        SwingDuration_SD = nanstd(SwingDuration);
        SwingDurationPercentage_SD = nanstd(SwingDurationPercentage);
        % Stance
        StanceDuration(Index_ExcSwingDuration) = [];
        StanceDurationPercentage(Index_ExcSwingDuration) = [];
        StanceDuration_CoeffVariation = nanstd(StanceDuration)/nanmean(StanceDuration);
        StanceDurationPercentage_CoeffVariation = nanstd(StanceDurationPercentage)/nanmean(StanceDurationPercentage);
        StanceDuration_SD = nanstd(StanceDuration);
        StanceDurationPercentage_SD = nanstd(StanceDurationPercentage);
    else %not enough steps to define average values
        % Swing
        SwingDuration_CoeffVariation = nan;
        SwingDurationPercentage_CoeffVariation = nan;
        SwingDuration_SD = nan;
        SwingDurationPercentage_SD = nan;
        % Stance
        StanceDuration_CoeffVariation = nan;
        StanceDurationPercentage_CoeffVariation = nan;
        StanceDuration_SD = nan;
        StanceDurationPercentage_SD = nan;        
    end %if StepsNumber
end %Settings
%end

% Outputs
% Swing
WB_SDMO.SwingDuration_CoeffVariation = nanmean(SwingDuration_CoeffVariation);
WB_SDMO.SwingDurationPercentage_CoeffVariation = nanmean(SwingDurationPercentage_CoeffVariation);
WB_SDMO.SwingDuration_SD = nanmean(SwingDuration_SD);
WB_SDMO.SwingDurationPercentage_SD = nanmean(SwingDurationPercentage_SD);
% Stance
WB_SDMO.StanceDuration_CoeffVariation = nanmean(StanceDuration_CoeffVariation);
WB_SDMO.StanceDurationPercentage_CoeffVariation = nanmean(StanceDurationPercentage_CoeffVariation);
WB_SDMO.StanceDuration_SD = nanmean(StanceDuration_SD);
WB_SDMO.StanceDurationPercentage_SD = nanmean(StanceDurationPercentage_SD);
end

