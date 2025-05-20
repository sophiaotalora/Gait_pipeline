function [WB_SDMO] = extract_PhaseFeatures_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)
%[WB_SDMO] = extract_PhaseFeatures_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings)

%% Calculation of phase plot features - ACCELERATION
% * Ref 1) According to: Michael Dunne-Willows, Professor Paul Watson, Dr Jian Shi, Dr Silvia Del Din
% "A Novel Parameterisation of Phase Plots for Monitoring of Parkinson’s Disease"
% https://eprint.ncl.ac.uk/258618
% "41st IEEE Engineering in Medicine and Biology Society conference", 2019
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
% WB_SDMO(i).PhaseFeatures_Acc = Phase plot features obtained from accelerometry data, for each signal-axis (the original script was based primarly on VT Acceleration)

%% Remarks
% *** Comments located on the right side
% *** Author: Encarna Micó Amigo, based on UNEW algorithm (BAM department, original author: Michael Dunne-Willows)
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

WB_Number = size(WB_WBD,2);                                                 % number of walking bouts for loop
for WB_index = 1:WB_Number                                                  % loop over each WB
    
    % Info about straight line walks within WB
    StraightWalk_Number = size(WB_WBD(WB_index).straightWalk,2);            % number/amount of straight walk episodes
    straightWalk_start = cell2mat({WB_WBD(WB_index).straightWalk.start});   % start of straight WB
    straightWalk_end = cell2mat({WB_WBD(WB_index).straightWalk.end});       % end of straight WB
    NameChannels = {'VT','ML','AP','Res'};                                  % "Res" instead of "Combined" as requested by the phase features function
    
    % Initialization
    clear PhaseFeatures
    VariablesPhaseNames = {'f_ell_ecc','f_ell_ang_asy','f_ell_area','f_ell_area_asy_prop','f_ell_gof','f_ell_minorSD','f_ell_majorSD','s_ell_ecc_asy','s_ell_ang_asy','s_ell_area_asy_prop','s_ell_gof','l_ell_ecc_asy','l_ell_ang_asy','l_ell_area_asy_prop','l_ell_gof','cor','gof_SD'};
    for iChannel = 1:size(NameChannels,2)
        for iPhase = 1:size(VariablesPhaseNames,1)
            PhaseFeatures.(NameChannels{iChannel}).(VariablesPhaseNames{iPhase}) = [];
        end
    end
    for StraightWalk_Index = 1:StraightWalk_Number                          % loop over straight walks within each WB
        
        % Definition of signal: start-end
        Acc = zeros(size(imu.acc_V(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index)),1),4);
        Acc(:,1) = imu.acc_V(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
        Acc(:,2) = imu.acc_ML(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
        Acc(:,3) = imu.acc_AP(straightWalk_start(StraightWalk_Index):straightWalk_end(StraightWalk_Index));
        if mean(Acc(:,1)) < 0                                               % signal negative = VT acceleration input is downwards
            Acc(:,1) = -Acc(:,1);                                           % we invert the signal (to get an upwards direction). If it would be directly provided as upwards, there would not be need to invert it
            Acc(:,3) = -Acc(:,3);                                           % we invert the signal (to get an forwards direction). If it would be directly provided as forwards, there would not be need to invert it
        end %although there is no need no invert it, since the Phase analysis function does it itself
        
        % Filter the signal (do not remove the trend in this case, since the Phase function requires to know the relative positon of the axes)
        Acc(:,4) = sqrt(sum(Acc(:,1:3)'.^2)');                              % resultant or combined signal, non-detrended because the function does it 
        Acc_LPFilt = WintFilt_low(Acc,Settings.Fc_LowPassFilter, fs);       % low-pass filtered
        
        % Initial temporal features calculation
        if ~isfield(WB_WBD,'TempEvents')                                     % not provided as an input (FAU & EPFL work)
            AccVT = Acc(:,1);
            if Settings.StepDetectMcCameley_EMA
                TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,fs,PlotSettings);
            else
                TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_BAM(Settings,AccVT,fs,PlotSettings);
            end
            IC = TemporalFeatures.InitialContactSample;
            StepDuration = TemporalFeatures.StepDuration;
            StepsNumber = size(IC,2);
        else                                                                % provided as an input
            IC = WB_WBD(WB_index).straightWalk.TempEvents;                  % [sample number]
            StepDuration = diff(IC)'./fs;                                   % [s]
            StepsNumber = WB_WBD(WB_index).straightWalk.stepsNumber;        % size(IC,2);
        end
        
        % Calculation of phase features calculation - Acc
        if  StepsNumber > Settings.MinNumStrides_Complex
            TimeVect = ones(size(Acc_LPFilt,1),1);
            AccSigInput = [TimeVect,Acc_LPFilt]; %signal should not be detrended, so the function for the analysis of phase plots can identify the VT axis by itself
            for iChannel = 1:size(Acc_LPFilt,2)
                if PlotSettings.PhaseFeatures
                    PhaseFeaturesPlotDelay = 2;
                else
                    PhaseFeaturesPlotDelay = 0;
                end
                PhaseFeatures(StraightWalk_Index).(NameChannels{iChannel}) = Phase_analysis_Acc(AccSigInput, NameChannelsPhase{iChannel}, fs, Settings.MinNumStrides_Complex, PlotSettings.PhaseFeatures, PhaseFeaturesPlotDelay);
            end
        else %not enough steps to define average values
            for iChannel = 1:size(Acc_LPFilt,2)
                VariablesPhaseNames = {'f_ell_ecc','f_ell_ang_asy','f_ell_area','f_ell_area_asy_prop','f_ell_gof','f_ell_minorSD','f_ell_majorSD','s_ell_ecc_asy','s_ell_ang_asy','s_ell_area_asy_prop','s_ell_gof','l_ell_ecc_asy','l_ell_ang_asy','l_ell_area_asy_prop','l_ell_gof','cor','gof_SD'};
                for iPhase = 1:size(VariablesPhaseNames,1)
                    PhaseFeatures(StraightWalk_Index).(NameChannels{iChannel}).(VariablesPhaseNames{iPhase}) = nan;
                end
            end
        end
    end
    
    % Outputs
    VariablesPhaseNames = fieldnames(PhaseFeatures(StraightWalk_Index).(NameChannels{iChannel})); %it will take the variable names (list) from the last straight walk bout and the last channel
    for iChannel = 1:size(Acc_LPFilt,2)
        for iPhase = 1:size(VariablesPhaseNames,1)
            DataStraightWalk = [];
            for iStraightWalks = 1:StraightWalk_Number %size(SpectralFeatures,2) % we need to accumulate all variables for each straight walk within the WB, and then calculate the mean value
                DataStraightWalk_ = PhaseFeatures(iStraightWalks).(NameChannels{iChannel}).(VariablesPhaseNames{iPhase}');
                DataStraightWalk = [DataStraightWalk; DataStraightWalk_];
            end
            WB_SDMO(WB_index).PhaseFeatures_Acc.(NameChannels{iChannel}).(VariablesPhaseNames{iPhase}) = nanmean(DataStraightWalk);
            clear DataStraightWalk_
        end
    end
end

