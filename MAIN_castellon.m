clear all;
%close all;
clc;

% Make sure your current working directory is the directory of this file
currentwd = pwd;
addpath(genpath('.'))
% Add folder with algorithm code to path
librarydir = absPath('C:\Users\sophi\Heriot-Watt University\Micó Amigo, Encarna - Digital Mobility Technology\Data\mobilise-d\algorithms');
addpath(genpath(librarydir));
%
% Set general pipeline parameters

% Define environment variables
SENSOR_UNIT     = 'SU'; % SU: sensor unit
SENSOR_POSITION = 'LowerBack'; %'LowerBack';  %'Wrist';

% Those variable are used to define the SU and Position names in the
% results output files. E.g., set them to 'SU', 'LowerBack', if you want to
% have the output struct fields always set to those values
% Default: Use the same as SENSOR_UNIT and SENSOR_POSITION
SENSOR_UNIT_OUTPUT_NAME = SENSOR_UNIT;
SENSOR_POSITION_OUTPUT_NAME = SENSOR_POSITION;

% Reference information
STANDARD_UNIT = 'IMU'; % Choose reference system: INDIP, Stereophoto, Walkway, IMU, Gaitrite,SU_LowerShanks
STANDARD_BOUT = 'MicroWB'; % Choose reference bout type: MicroWB, ContinuousWalkingPeriod, Pass,

% Plot intermediate results?
BLOCK_PLOT = 'false';

setenv('SENSOR', SENSOR_UNIT);
setenv('SENSOR_POSITION',SENSOR_POSITION);
setenv('SENSOR_UNIT_OUTPUT_NAME', SENSOR_UNIT_OUTPUT_NAME);
setenv('SENSOR_POSITION_OUTPUT_NAME',SENSOR_POSITION_OUTPUT_NAME);
setenv('STANDARD',STANDARD_UNIT);
setenv('BOUT',STANDARD_BOUT);
setenv('PLOT', BLOCK_PLOT);
% [0] PREPARATION & ACTIVATION & SETTINGS
Fs=100;
% 0.4 Settings for features extraction
for Close_SetSettings = 1                                               % close settings
    Settings = struct;                                                  % preparation of a structure for general settings
    % General
    Settings.UseWalkingIdentification =  1;                             % use "walking activities" identified or the complete signal
    % Specific
    Settings.SigLevel = 0.05;                                           % alpha value for ttest, scalar between 0 and 1. Default = 0.05
    Settings.Fs = 100; %[Hz]
    Settings.Fc_LowPassFilter = 20;                                     % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter_AutoCorr = 0.8;                          % cut-off frequency used for high-pass filters for autocorrelation unbiased: : frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.N_Harm = 20;                                               % number of harmonics used for harmonic ratio, index of harmonicity and phase fluctuation
    Settings.MaxRangeFrequencyPSD = 10;                                 % maximal range of frequencies to find frequency of the peak power spectral density (PSD)
    Settings.MinRangeFrequencyPSD = 0.3;                                % miniaml range of frequencies to find frequency of the peak power spectral density (PSD)
    Settings.JunkSizes = [3 5 7 10]; %[s]                               % durations of windows for spectra features extraction
    % Settings for Spatiotemporal extraction
    %Settings.MinNumEvents = eval(getenv(‘MIN_NUM_EVENTS’));
    Settings.MinNumEvents = 3;                  %CONSENSUS!             % minimal number of steps to consider the step-to-stpe gait features
    Settings.PercOfDominantFreq = 0.5;                                  % percentage of the dominant frequency considered as a threshold (minimal distance between peaks)
    Settings.GravityConstant = 9.81; % [m/s2]     
    %mirar filtros
    Settings.Fc_LowPassFilter = 20*Fs/100;                              % cut-off frequency used for low-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter = 0.1*Fs/100;                            % cut-off frequency used for high-pass filters: frequencies scale according to the sampling rate, so it must be scaled by 1.28, since most of the Fc are defined in the literature for a Fs of 100 Hz
    Settings.Fc_HighPassFilter_AutoCorr = 0.8*Fs/100;                   % cut-off frequency used for high-pass filters, particularly tailored for the autocorrelation signal of the 1st derivative 
    Settings.StepDetectMcCameley_EMA = 0;                               % use McCamely algorithm with Encarna modifications or with BAM original 
    % Settings for step length calculation
    Settings.L5_height = 1.05; % [m]                                    % height of the sensor: from the floor to L5 (particular to the subject)
    %REVISAR VALORES EN PAPER ACTUAL
    Settings.CorrectionFactor = 1.25;                                   % correction factor defined by Zijlstra, based no the iverted pendulum model
    Settings.MaxRangeStepLength = 0.85;         %CONSENSUS!             % maximal range for step length
    Settings.MinRangeStepLength = 0.23;         %CONSENSUS!             % minimal range for step length
    Settings.MaxRangeStepDuration = 1.25;       %CONSENSUS!             % maximal range for step duration
    Settings.MinRangeStepDuration = 0.25;       %CONSENSUS!             % minimal range for step duration
    Settings.MaxRangeSwingDuration = 0.85;      %CONSENSUS!             % maximal range for swing duration
    Settings.MinRangeSwingDuration = 0.23;      %CONSENSUS!             % minimal range for swing duration
    Settings.StrideDurationRange = [0.4 4.0];                           % Range to search for stride duration (seconds)
    % Settings for frequency analysis
    Settings.N_Harm = 20;                                               % number of harmonics used for harmonic ratio, index of harmonicity and phase fluctuation
    Settings.MaxRangeFrequencyPSD = 10*Fs/100;                          % maximal range of frequencies to find frequency of the peak power spectral density (PSD)
    Settings.MinRangeFrequencyPSD = 0.3*Fs/100;                         % miniaml range of frequencies to find frequency of the peak power spectral density (PSD)
    % Settings for phase/complex measures
    Settings.MinNumStrides_Complex = 7;%floor(Settings.MinNumEvents/2); % minimal number of strides in the signal
    Settings.HanningSize = 4;  
    Settings.GyrUnits = 1;                                              %units of gyroscope. 1= deg/s, 0 = rad/s
    Settings.Given_data = 1;                                            % 1 = Given data (excel), 0 = use elgohary algorithm 
end
for Close_PlotSettings = 1
    %PlotSettings.PhaseFeatures = 0;                                    % activate plots to help understangind the development of the phase analysis algorithms
    PlotSettings.SpectralFatures = 0;                                   % activate plots to help understangind the calculation of power spectral density of signals      
    PlotSettings.StepDetectionDevelopment = 0;                          % activate plots to help understangind the development of the step detection algorithm
    PlotSettings.StepLengthDevelopment = 0;                             % activate plots to help understangind the development of the step length calculation algorithm
    PlotSettings.PhaseFeatures = 0;                                     % activate plots to help understangind the development of the phase analysis algorithms
end
 
clear Close_Prepare_Settings Close_Initialize_Variables_Plots Close_Define_Directory Close_ActivateAnalyses ...
    Directory_Catalyst Close_SetSettings Close_PlotSettings

%%
% Create data.mat
DirectoryPath_Data = ('C:\Users\sophi\Heriot-Watt University\Micó Amigo, Encarna - Digital Mobility Technology\Data\Castellon_data');
DirectoryPath_Functions= ('C:\Users\sophi\Heriot-Watt University\Micó Amigo, Encarna - Digital Mobility Technology\Data\Castellon_data\Functional Functions');
cd(DirectoryPath_Data)
fileName= 'data.mat';
fullSavePath = fullfile(DirectoryPath_Data, fileName);
indir = DirectoryPath_Data; %cambiar sophi por so3003
outdir = indir;
load(fullfile(outdir, 'data.mat'), 'Data');
[Data]= add_imufield(Data);

%% Gait Sequence detection GSD
close all;
cd(fullfile(librarydir, 'GSDA'));
driver(indir, outdir);
cd(currentwd);
%
load(fullfile(outdir, 'GSDA_Output.mat'), 'GSD_Output');

%% Turns
% Turns
close all;
%imu= data.TimeMeasure1.Test1.Trial1.SU;
cd(fullfile(librarydir, 'TD_ElGohary'));
driver(indir, outdir);
cd(currentwd);
%
load(fullfile(outdir, 'TD_ElGohary_SU.mat'), 'TD_Output');
%% Divide signal straight and turns
[Data] = divide_signal_c(Data, TD_Output);
%% extract features
[Data] = extract_features_c(Fs, Data, Settings, PlotSettings);
%% Make excel with mean value
load(fullfile(outdir, 'Data_final.mat'), 'Data');
%%
States = fieldnamesr(Data, 2);
States = unique(regexprep(States, '.*\.', ''));
SubjectNames = fieldnames(Data);
Trials= fieldnamesr(Data, 3);
Trials = unique(regexprep(Trials, '.*\.', ''));
for iSubject = 1:length(SubjectNames)
    for iState = 1:length(States)
        for iTrial = 1:length(Trials)
            Features= fieldnames(Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO);
            Features_size= size(Features,1);
            for iFeature = 1:Features_size
                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.(Features{iFeature}) = nonzeros(Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.(Features{iFeature}))';
                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.(Features{iFeature}) = nanmean(Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.(Features{iFeature})); 
                if Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.(Features{iFeature})== 25
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO = rmfield(Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO, Features{iFeature});
                end
            end
            %erase zeros
            %WB_SDMO=  Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO;
            
        end
    end
end





