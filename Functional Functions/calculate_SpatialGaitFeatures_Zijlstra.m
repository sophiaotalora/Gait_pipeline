function SpatialFeatures = calculate_SpatialGaitFeatures_Zijlstra(AccVT,Fs,IC,StepDuration,Settings,PlotSettings)
%function SpatialFeatures = calculate_SpatialGaitFeatures_Zijlstra(AccVT,Fs,IC,StepDuration,Settings,PlotSettings)
%
% OUTPUTS
% * TemporalFeatures.
%       - StepLength = displacement covered for every step;
%       - StepVelocity = velocity of each step;
%
% AUTHOR
% Encarna Micó Amigo (based on algorithm developed at Newcastle University - BAM department, Silvia Del Din)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
% Objective: Development of an algorithm to calculate walking speed related features: walking speed and stpe length
%
% GENERAL DESCRIPTION
% Calculate stpe length and step velocity from segmented accelerometry 
% According to: Zijlstra, Wiebren.
% "Assessment of spatio-temporal parameters during unconstrained walking."
% European journal of applied physiology 92,(2004): 39-44"
%
% INPUT
% * AccVT  = triaxial acceleration of the defined gait episode (vertical axis), downwards (negative) / upwards (positive)
% * Fs = sampling frequency
% * IC = initial contact or heel strike
% * StepDuration = duration of step cycles
% * Settings. 
%        - L5_height: minimal number of events required to extract temporal gait features
%        - GravityConstant: gravity acceleration
%        - Fc_HighPassFilter: cutting frequency for a high-pass filter
%        - CorrectionFactor: correction factor required by the model Zijlstra

% * PlotSettings
%        - StepLengthDevelopment: 1 (activate plotting), 0 (deactivate)
%       
% HISTORY
% 2019/3rd/September functionized - Encarna Micó Amigo
%
% REMARKS
% * Comments on the right side


%% [0] DEFINITION OF VARIABLES & PRE-PROCESSING

% 0.1 Definition
SensorHeight = Settings.L5_height;                                          % distance from the floor to the L5 (lower-back) from the participant
GravityConstant = Settings.GravityConstant;                                 % Units: m/second2
Fc_LowPassFilter = Settings.Fc_LowPassFilter;                               % cut off frequency used for the low-pass filter
Fc_HighPassFilter = Settings.Fc_HighPassFilter;                             % cut off frequency used for the high-pass filter
% 0.2 Pre-Processing
if mean(AccVT) < 0                                                          % signal negative = VT acceleration input is downwards
    AccVT = -AccVT;                                                         % the algorithm requires a positive direction of VT Acc, thus, we ivert the signal (as upwards). If it would be directly provided as upwards, there would not be need to invert it                                                                                        
end
AccVT_Detr = detrend(AccVT);                                                % remove the trend/offset of gravity
AccVT_Detr_Filt = WintFilt_low(AccVT_Detr-mean(AccVT_Detr),Fc_LowPassFilter,Fs); %low-pass filter the signal


%% [1] DOUBLE INTEGRATION OF ACCELERATION ---> DISTANCE
 
% 1.1 Double integration
VelocityVT = cumtrapz(AccVT_Detr_Filt).*1/Fs;
VelocityVT_Detr = detrend(VelocityVT); 
PositionVT = cumtrapz(VelocityVT_Detr).*1/Fs;
PositionVT_Detr = detrend(PositionVT);  
PositionVT_Detr_Filt = WintFilt_high(PositionVT_Detr,Fc_HighPassFilter,Fs); % filter 4th order buttherworth Filter

% 1.2 Position per step
PositionAtInitialContacts = PositionVT_Detr_Filt(IC);
for iIC = 1:size(IC,2)-1
    PositionPerStep = PositionVT_Detr_Filt(IC(iIC):IC(iIC+1));
    HeightChangePerStep(iIC,1) = abs(max(PositionPerStep) - min(PositionPerStep));    
end
%try with one step
% PositionAtInitialContacts = PositionVT_Detr_Filt(IC);
% if size(IC,2)==2
%     for iIC = 1:size(IC,2)-1
%         PositionPerStep = PositionVT_Detr_Filt(IC(iIC):IC(iIC+1));
%         HeightChangePerStep(iIC,1) = abs(max(PositionPerStep) - min(PositionPerStep));    
%     end
% else
%     for iIC = 1:size(IC,2)
%         PositionPerStep = PositionVT_Detr_Filt(IC(iIC));
%         HeightChangePerStep(iIC,1) = PositionPerStep;    
%     end
% end

% 1.3 Inverted Pendulum Zijlstra Model
StepLength_Incorrected = 2.*sqrt(abs(2.*SensorHeight.*HeightChangePerStep - HeightChangePerStep.^2));
StepLength = StepLength_Incorrected.*Settings.CorrectionFactor;             % based on the paper, the obtained outcome must be scaled by a correction factor of 1.25
StepVelocity = StepLength./StepDuration;


%% [2] PLOTS
if PlotSettings.StepLengthDevelopment
    figure, hold on, plot(VelocityVT_Detr,'m'), plot(PositionVT_Detr_Filt,'k')
    %plot(AccVT_Detr_Filt,'b'),
    ax=axis;
    line([IC(:) IC(:)],[ax(3) ax(4)],'Color' ,'c','LineWidth',0.9);
end


%% [3] OUTPUTS

SpatialFeatures.StepLength = StepLength;
SpatialFeatures.StepVelocity = StepVelocity;
 

