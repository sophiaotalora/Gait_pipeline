function [WB_SDMO] = extract_AngularAccMax_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_AngularAccMax_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of Root-mean-square of signals (pure) - 1st DERIVATIVE OF GYROSCOPE DATA / ANGULAR ACCELERATION
% * Ref 1) According to: Luca Palmerini, Sabato Mellone, Guido Avanzolini, Franco Valzania, and Lorenzo Chiari.
% "Quantification of Motor Impairment in Parkinson’s Disease Using an Instrumented Timed Up and Go Test"
% IEEE TRANSACTIONS ON NEURAL SYSTEMS AND REHABILITATION ENGINEERING, VOL. 21, NO. 4, JULY 2013
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
% WB_SDMO(i).AngAccMax_SigComplete = Max of angular acceleration (1st derivative of angular velocity) of complete signal (walking bout in the given imu data).
% WB_SDMO(i).AngAccMax_AvPerSteps = average Max of angular acceleration (1st derivative of angular velocity) of all steps within the WB (walking bout in the given imu data).
% WB_SDMO(i).AngAccMax_AvPerStrides = average Max of angular acceleration (1st derivative of angular velocity) of all stride within the WB (walking bout in the given imu data).

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

WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
for WB_index = 1:WB_Number                                                  % loop over each WB
    
    % Info about straight line walks within WB
    StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);            % number/amount of straight walk episodes
    straightWalk_start = cell2mat({WB_WBD(WB_index).straightWalk.start});   % start of straight WB
    straightWalk_end = cell2mat({WB_WBD(WB_index).straightWalk.end});       % end of straight WB
    NameChannels = {'Yaw','Pitch','Roll','Combined'};
    
    % Initialization
    clear AngAccMax_SigComplete
    clear AngAccMax_AvPerSteps
    clear AngAccMax_AvPerStrides
    for iChannel = 1:size(NameChannels,2)
        AngAccMax_SigComplete.(NameChannels{iChannel}) = [];
        AngAccMax_AvPerSteps.(NameChannels{iChannel}) = [];
        AngAccMax_AvPerStrides.(NameChannels{iChannel}) = [];
    end
    for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
        
        % Definition of signal: start-end
        Gyr = zeros(size(imu.gyro_yaw(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index)),1),4);
        Gyr(:,1) = imu.gyro_yaw(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
        Gyr(:,2) = imu.gyro_pitch(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
        Gyr(:,3) = imu.gyro_roll(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
        
        % Detrend & filter signals
        Gyr_Detr = detrend(Gyr);                                            % remove the trend
        Gyr_Detr_LPFilt = WintFilt_low(Gyr_Detr-mean(Gyr_Detr),Settings.Fc_LowPassFilter, fs); %low-pass filter the signal
        Gyr_Detr_LPFilt(:,4) = sqrt(sum(Gyr_Detr_LPFilt(:,1:3)'.^2)');      % resultant or combined signal
        
        % Initial temporal features calculation
        if ~isfield(WB_WBD,'TempEvents')                                    % not provided as an input (FAU & EPFL work)
            AccVT = imu.acc_V(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
            if Settings.StepDetectMcCameley_EMA
                TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
            else
                TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
            end
            IC = TemporalFeatures.InitialContactSample;
            StepsNumber = size(IC,2);
        else                                                                % provided as an input
            IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
            StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);
        end
        
        % Calculation of Max of Angular Acceleration - Gyro
        if StepsNumber > Settings.MinNumEvents
            for iChannel = 1:size(Gyr_Detr_LPFilt,2)
                AngAccMax_SigComplete(StraightWalk_Index).(NameChannels{iChannel}) = max(diff(Gyr_Detr_LPFilt(:,iChannel))*fs);% Angular acceleration is the 1st derivative of angular velocity, over it, we calcuate max
                for iIC = 1:size(IC,2)-1
                    AngAccMax_PerSteps(iIC,iChannel) = max(diff(Gyr_Detr_LPFilt(IC(iIC):IC(iIC+1),iChannel))*fs);
                    if iIC < size(IC,2)-1
                        AngAccMax_PerStrides(iIC,iChannel) = max(diff(Gyr_Detr_LPFilt(IC(iIC):IC(iIC+2),iChannel))*fs);
                    end
                end
                AngAccMax_AvPerSteps(StraightWalk_Index).(NameChannels{iChannel}) = nanmean(AngAccMax_PerSteps(:,iChannel));
                AngAccMax_AvPerStrides(StraightWalk_Index).(NameChannels{iChannel}) = nanmean(AngAccMax_PerStrides(:,iChannel));
            end
        else
            for iChannel = 1:size(Gyr_Detr_LPFilt,2)
                AngAccMax_SigComplete(StraightWalk_Index).(NameChannels{iChannel}) = nan;
                AngAccMax_AvPerSteps(StraightWalk_Index).(NameChannels{iChannel}) = nan;
                AngAccMax_AvPerStrides(StraightWalk_Index).(NameChannels{iChannel}) = nan;
            end
        end
    end
    
    % Outputs
    for iChannel = 1:size(Gyr_Detr_LPFilt,2)
        WB_SDMO(WB_index).AngAccMax_SigComplete.(NameChannels{iChannel}) = nanmean(extractfield(AngAccMax_SigComplete,(NameChannels{iChannel})));
        WB_SDMO(WB_index).AngAccMax_AvPerSteps.(NameChannels{iChannel}) = nanmean(extractfield(AngAccMax_AvPerSteps,(NameChannels{iChannel})));
        WB_SDMO(WB_index).AngAccMax_AvPerStrides.(NameChannels{iChannel}) = nanmean(extractfield(AngAccMax_AvPerStrides,(NameChannels{iChannel})));
    end    
end

