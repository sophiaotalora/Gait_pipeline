function [WB_SDMO] = extract_JerkMeanLogRatio_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_JerkMeanLogRatio_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Root-mean-square of signals (pure) - 1st DERIVATIVE OF ACCELEROMETRY
% * Ref 1) According to: Brodie MA, Menz HB, Lord SR.
% "Age-associated changes in head jerk while walking reveal altered dynamic stability in older people"
% Exp Brain Res. 2014 Jan;232(1):51-60. doi: 10.1007/s00221-013-3719-6. Epub 2013 Oct 5
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

%% Output
% WB_SDMO(i).JerkLogRatio_SigComplete = normalized Ratio with logarithmic transformation of jerk (1st derivative of acceleration) of complete signal (walking bout in the given imu data),i.e. jerk on each component / jerk on resultant
% WB_SDMO(i).JerkLogRatio_AvPerSteps = average of normalized ratio with logarithmic transformation of jerk (1st derivative of acceleration) of all steps within the WB (walking bout in the given imu data), i.e. jerk on each component / jerk on resultant
% WB_SDMO(i).JerkLogRatio_AvPerStrides = average of normalized ratio with logarithmic transformation of jerk (1st derivative of acceleration) of all stride within the WB (walking bout in the given imu data), i.e. jerk on each component / jerk on resultant

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Mic� Amigo, based on UNEW algorithm (BAM department, original author: Chris Buckley)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
     
%% History
% 2019/6th/September functionized - Encarna Mic� Amigo


%% SETTINGS

if nargin < 5                                                               % if the settings are not predefined:
    % Settings for Spatiotemporal extraction
    Settings.MinNumEvents = 3;                  %CONSENSUS!                 % minimal number of steps to consider the step-to-stpe gait features
    Settings.GravityConstant = 9.81; % [m/s2]
    Settings.Fc_LowPassFilter = 20*fs/100;                                  % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter = 0.1*fs/100;                                % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter_AutoCorr = 0.8*fs/100;                       % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.StepDetectMcCameley_EMA = 1;                                   % use McCamely algorithm with Encarna modifications (1) or with BAM original (0)
end

%% LOOP OVER WALKING BOUTS

%WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
%for WB_index = 1:WB_Number                                                  % loop over each WB
    
% Info about straight line walks within WB
%StraightWalk_Number = size(WB_WBD.straightWalk,2);            % number/amount of straight walk episodes
straightWalk_start = cell2mat({WB_WBD.start});   % start of straight WB
straightWalk_end = cell2mat({WB_WBD.end});       % end of straight WB
NameChannels = {'VT','ML','AP','Combined'};

% Initialization
clear JerkLogRatio_SigComplete
clear JerkLogRatio_AvPerSteps
clear JerkLogRatio_AvPerStrides
for iChannel = 2:3                                                      % only ML and AP %%why??
    JerkLogRatio_SigComplete.(NameChannels{iChannel}) = [];
    JerkLogRatio_AvPerSteps.(NameChannels{iChannel}) = [];
    JerkLogRatio_AvPerStrides.(NameChannels{iChannel}) = [];
end
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
    
% Definition of signal: start-end
Acc = zeros(size(imu.acc_V(straightWalk_start:straightWalk_end),1),3);
Acc(:,1) = imu.acc_V(straightWalk_start:straightWalk_end);
Acc(:,2) = imu.acc_ML(straightWalk_start:straightWalk_end);
Acc(:,3) = imu.acc_AP(straightWalk_start:straightWalk_end);
if mean(Acc(:,1)) < 0                                               % signal negative = VT acceleration input is downwards
    Acc(:,1) = -Acc(:,1);                                           % we invert the signal (to get an upwards direction). If it would be directly provided as upwards, there would not be need to invert it
    Acc(:,3) = -Acc(:,3);                                           % we invert the signal (to get an forwards direction). If it would be directly provided as forwards, there would not be need to invert it
end % although there is no need no invert it

% Detrend & filter signals
Acc_Detr = detrend(Acc);                                            % remove the trend/offset of gravity
Acc_Detr_LPFilt = WintFilt_low(Acc_Detr-mean(Acc_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal

% Initial temporal features calculation
%         if ~isfield(WB_WBD,'TempEvents')                                    % not provided as an input (FAU & EPFL work)
AccVT = imu.acc_V(straightWalk_start:straightWalk_end);
if Settings.StepDetectMcCameley_EMA
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
else
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
end
IC = TemporalFeatures.InitialContactSample;
StepsNumber = size(IC,2);
%         else                                                                % provided as an input
%             IC = WB_WBD.straightWalk.TempEvents;                  % [sample number]
%             StepsNumber = WB_WBD.straightWalk.stepsNumber;        % size(IC,2);
%         end

% Calculation of Mean Normalized Ratio of Jerk - Acc
if StepsNumber > Settings.MinNumEvents
    for iChannel = 2:3                                              % ML and AP
        JerkLogRatio_SigComplete.(NameChannels{iChannel}) = mean(10*log10((diff(Acc_Detr_LPFilt(:,iChannel))*fs)./(diff(Acc_Detr_LPFilt(:,1))*fs)));% Jerk is the 1st derivative of acceleration, then, we calcuate the normalized (with log) ratio: jerk on each component / jerk on resultant
        for iIC = 1:size(IC,2)-1
            JerkLogRatio_PerSteps(iIC,iChannel) = mean(10*log10((diff(Acc_Detr_LPFilt(IC(iIC):IC(iIC+1),iChannel))*fs)./(diff(Acc_Detr_LPFilt(IC(iIC):IC(iIC+1),1))*fs)));
            if iIC < size(IC,2)-1
                JerkLogRatio_PerStrides(iIC,iChannel) = mean(10*log10((diff(Acc_Detr_LPFilt(IC(iIC):IC(iIC+2),iChannel))*fs)./(diff(Acc_Detr_LPFilt(IC(iIC):IC(iIC+2),1))*fs)));
            end
        end
        JerkLogRatio_AvPerSteps.(NameChannels{iChannel}) = nanmean(JerkLogRatio_PerSteps(:,iChannel));
        JerkLogRatio_AvPerStrides.(NameChannels{iChannel}) = nanmean(JerkLogRatio_PerStrides(:,iChannel));
    end
else
    for iChannel = 2:3
        JerkLogRatio_SigComplete.(NameChannels{iChannel}) = nan;
        JerkLogRatio_AvPerSteps.(NameChannels{iChannel}) = nan;
        JerkLogRatio_AvPerStrides.(NameChannels{iChannel}) = nan;
    end
end
%end

% Outputs
for iChannel = 2:3
    WB_SDMO.JerkMeanLogRatio_SigComplete.(NameChannels{iChannel}) = nanmean(extractfield(JerkLogRatio_SigComplete,(NameChannels{iChannel})));
    %WB_SDMO.JerkMeanLogRatio_AvPerSteps.(NameChannels{iChannel}) = nanmean(extractfield(JerkLogRatio_AvPerSteps,(NameChannels{iChannel})));
    %WB_SDMO.JerkMeanLogRatio_AvPerStrides.(NameChannels{iChannel}) = nanmean(extractfield(JerkLogRatio_AvPerStrides,(NameChannels{iChannel})));
end    
%end

