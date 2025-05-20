function [WB_SDMO] = extract_RMSRatio_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_RMSRatio_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Root-mean-square of signals (ratio)  - GYROSCOPE / ANGULAR VELOCITY
% * Ref 1) According to: Sekine M, Tamura T, Yoshida M, Suda Y, Kimura Y, Miyoshi H, Kijima Y, Higashi Y, Fujimoto T..
% "A gait abnormality measure based on root mean square of trunk acceleration."
% J Neuroeng Rehabil,(2013): https://www.ncbi.nlm.nih.gov/pubmed/24370075
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
% WB_SDMO(i).RMSratio_SigComplete = RMS ratio of complete signal (walking bout in the given imu data), i.e. RMS of each component / RMS of resultant
% WB_SDMO(i).RMSratio_AvPerSteps = average RMS ratio of all steps within the WB (walking bout in the given imu data), i.e. RMS of each component / RMS of resultant
% WB_SDMO(i).RMSratio_AvPerStrides = average RMS ratio of all stride within the WB (walking bout in the given imu data), i.e. RMS of each component / RMS of resultant

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Micó Amigo, Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
      
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
clear RMSratio_SigComplete
clear RMSratio_AvPerSteps
clear RMSratio_AvPerStrides
for iChannel = 1:size(NameChannels,2)-1
    RMSratio_SigComplete.(NameChannels{iChannel}) = [];
    RMSratio_AvPerSteps.(NameChannels{iChannel}) = [];
    RMSratio_AvPerStrides.(NameChannels{iChannel}) = [];
end
%for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
    
% Definition of signal: start-end
Gyr = zeros(size(imu.acc_V(straightWalk_start:straightWalk_end),1),4);
Gyr(:,1) = imu.gyro_yaw(straightWalk_start:straightWalk_end);
Gyr(:,2) = imu.gyro_pitch(straightWalk_start:straightWalk_end);
Gyr(:,3) = imu.gyro_roll(straightWalk_start:straightWalk_end);

% Detrend & filter signals
Gyr_Detr = detrend(Gyr);                                            % remove the trend
Gyr_Detr_LPFilt = WintFilt_low(Gyr_Detr-mean(Gyr_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal
Gyr_Detr_LPFilt(:,4) = sqrt(sum(Gyr_Detr_LPFilt(:,1:3)'.^2)');      % resultant or combined signal

% Initial temporal features calculation
% if ~isfield(WB_WBD,'TempEvents')                                    % not provided as an input (FAU & EPFL work)
AccVT = imu.acc_V(straightWalk_start:straightWalk_end);
if Settings.StepDetectMcCameley_EMA
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
else
    TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
end
IC = TemporalFeatures.InitialContactSample;
StepsNumber = size(IC,2);
% else                                                                % provided as an input
%     IC = WB_WBD.straightWalk.TempEvents;                  % [sample number]
%     StepsNumber = WB_WBD.straightWalk.stepsNumber;        % size(IC,2);
% end
        
% Calculation of RMS Ratio - Gyro
if StepsNumber > Settings.MinNumEvents
    for iChannel = 1:size(Gyr_Detr_LPFilt,2)-1
        RMSratio_SigComplete.(NameChannels{iChannel}) = rms(Gyr_Detr_LPFilt(:,iChannel))./rms(Gyr_Detr_LPFilt(:,4));
        for iIC = 1:size(IC,2)-1
            RMSratio_PerSteps(iIC,iChannel) = rms(Gyr_Detr_LPFilt(IC(iIC):IC(iIC+1),iChannel))./rms(Gyr_Detr_LPFilt(IC(iIC):IC(iIC+1),4));
            if iIC < size(IC,2)-1
                RMSratio_PerStrides(iIC,iChannel) = rms(Gyr_Detr_LPFilt(IC(iIC):IC(iIC+2),iChannel))./rms(Gyr_Detr_LPFilt(IC(iIC):IC(iIC+2),4));
            end
        end
        RMSratio_AvPerSteps.(NameChannels{iChannel}) = nanmean(RMSratio_PerSteps(:,iChannel));
        RMSratio_AvPerStrides.(NameChannels{iChannel}) = nanmean(RMSratio_PerStrides(:,iChannel));
    end
else
    for iChannel = 1:size(Gyr_Detr_LPFilt,2)-1
        RMSratio_SigComplete.(NameChannels{iChannel}) = nan;
        RMSratio_AvPerSteps.(NameChannels{iChannel}) = nan;
        RMSratio_AvPerStrides.(NameChannels{iChannel}) = nan;
    end
end
%end

% Outputs
for iChannel = 1:size(Gyr_Detr_LPFilt,2)-1
    WB_SDMO.RMSratio_SigComplete.(NameChannels{iChannel}) = nanmean(extractfield(RMSratio_SigComplete,(NameChannels{iChannel})));
    %WB_SDMO.RMSratio_AvPerSteps.(NameChannels{iChannel}) = nanmean(extractfield(RMSratio_AvPerSteps,(NameChannels{iChannel})));
    %WB_SDMO.RMSratio_AvPerStrides.(NameChannels{iChannel}) = nanmean(extractfield(RMSratio_AvPerStrides,(NameChannels{iChannel})));
end    
%end

