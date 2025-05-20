function [WB_SDMO] = extract_AllBouts_Variability(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_AllBoutsVariability(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of variability in bouts duration (for straight gait episodes within a single episode)
% * Ref 1) Authors: Rochester et al.
% Understanding the impact of deep brain stimulation on ambulatory activity in advanced Parkinson’s disease, 
% J Neurol, 2012
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
%       Settings.Fc_LowPassFilter = 20;                                     % cut-off frequency used for low-pass filters
%       Settings.Fc_HighPassFilter = 0.1;                                   % cut-off frequency used for high-pass filters
%       Settings for step length calculation
%       Settings.CorrectionFactor = 1.25;                                   % correction factor defined by Zijlstra, based no the iverted pendulum model
%       Settings.MaxRangeStepLength = 0.85;         %CONSENSUS!             % minimal range for step length
%       Settings.MinRangeStepLength = 0.23;         %CONSENSUS!             % maximal range for step length
% 6) Settings for Plots created (for the development/checking of algorithm)
%       PlotSettings.StepDetectionDevelopment = 0;
%       PlotSettings.StepLengthDevelopment = 0;

%% Output
% WB_SDMO(i).SD_WBDuration = standard deviation of walking bouts duration 
% WB_SDMO(i).CoeffVariation_WBDuration = coefficient of variation of walking bouts duration

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Micó Amigo, based on UNEW algorithm (BAM department, original author: Silvia Del Din)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
       
%% History
% 2019/6th/September functionized - Encarna Micó Amigo

 %% SETTINGS 

 if nargin < 5                                                              % if the settings are not predefined:
     %PlotSettings
     PlotSettings.StepDetectionDevelopment = 0;
     PlotSettings.StepLengthDevelopment = 0;
     %Settings for Spatiotemporal extraction
     Settings.MinNumEvents = 3;                  %CONSENSUS!                % minimal number of steps to consider the step-to-stpe gait features
     Settings.PercOfDominantFreq = 0.5;                                     % percentage of the dominant frequency considered as a threshold (minimal distance between peaks)
     Settings.L5_height = 1.05; % [m]                                       % height of the sensor: from the floor to L5
     Settings.GravityConstant = 9.81; % [m/s2]
     Settings.Fc_LowPassFilter = 20;                                        % cut-off frequency used for low-pass filters
     Settings.Fc_HighPassFilter = 0.1;                                      % cut-off frequency used for high-pass filters
     Settings.Fc_HighPassFilter_AutoCorr = 0.8;                             % cut-off frequency used for high-pass filters for autocorrelation unbiased
     %Settings for step length calculation
     Settings.CorrectionFactor = 1.25;                                      % correction factor defined by Zijlstra, based no the iverted pendulum model
     Settings.MaxRangeStepLength = 0.85;         %CONSENSUS!                % minimal range for step length
     Settings.MinRangeStepLength = 0.23;         %CONSENSUS!                % maximal range for step length
 end

%% LOOP OVER WALKING BOUTS

WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
for WB_index = 1:WB_Number                                                  % loop over each WB
    
    % Info about straight line walks within WB
    StriaghtWalk_Number = size(WB_WBD(WB_index).straightWalk,2);            % number/amount of straight walk episodes
    straightWalk_start = cell2mat({WB_WBD(WB_index).straightWalk.start});   % start of straight WB
    straightWalk_end = cell2mat({WB_WBD(WB_index).straightWalk.end});       % end of straight WB
    
    % Initialization
    WalkingBoutDuration = [];
    for StriaghtWalk_Index = 1:StriaghtWalk_Number                          % loop over straight walks within each WB
        
        % Definition of signal: start-end & Number of identified steps
        AccVT = imu.acc_V(straightWalk_start(StriaghtWalk_Index):straightWalk_end(StriaghtWalk_Index));
        AccML = imu.acc_ML(straightWalk_start(StriaghtWalk_Index):straightWalk_end(StriaghtWalk_Index));
        AccAP = imu.acc_AP(straightWalk_start(StriaghtWalk_Index):straightWalk_end(StriaghtWalk_Index));
                
        % Initial temporal features calculation
        if ~isfield(WB_WBD,'TempEvents')                                     % not provided as an input (FAU & EPFL work)
            TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
            %TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
            IC = TemporalFeatures.InitialContactSample;
            StepDuration = TemporalFeatures.StepDuration;
            StepsNumber = size(IC,2); 
        else                                                                % provided as an input 
            IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
            StepDuration = diff(IC)'./fs;                                   % [s]
            StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);            
        end
        
        % Walking Bout duration 
        if  StepsNumber > Settings.MinNumEvents
            WalkingBoutDuration(StriaghtWalk_Index) = (straightWalk_end(StriaghtWalk_Index)-straightWalk_start(StriaghtWalk_Index))./fs;  
        else %not enough steps to define average values
            WalkingBoutDuration(StriaghtWalk_Index) = nan;
        end
    end
    
    % Outputs: mean WBs duration, std and coefficient of variation
    WB_SDMO(WB_index).SD_WBDuration = nanstd(WalkingBoutDuration);
    WB_SDMO(WB_index).CoeffVariation_WBDuration = nanstd(WalkingBoutDuration)/nanmean(WalkingBoutDuration);
end

