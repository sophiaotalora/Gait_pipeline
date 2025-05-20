function TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,Fs,PlotSettings)
%function TemporalFeatures = calculate_TemporalGaitFeatures_McCamley_EMA(Settings,AccVT,Fs,PlotSettings)
%
% OUTPUTS
% * TemporalFeatures.
%       - InitialContactSample = samples from the signal identified as the initial contact or heel-strikes
%       - FinalContactsSample = samples from the signal identified as the final contact or toe-offs
%       - StepDuration = duration of step cycles, defined between initial contacts
%       - StrideDuration = duration of stride cycles, defined between two initial contacts (non consecutive, i.e. from the same foot)
%       - StanceDuration = duration between initial contact and final contact;
%       - SwingDuration = stride duration, except the stance duration;
%
% AUTHOR
% Encarna Micó Amigo (based on algorithm developed at Newcastle University - BAM department, original authors: Alan Godfrey & Silvia Del Din)
% Contact: Maria.Mico-Amigo@newcastle.ac.uk / encarna.mico@gmail.com
% Objective: Development of an algorithm for the detection of IC (initial contacts) and FC (Final Contact - toe off) and calculation of temporal gait features (related), such as: step duration, stride duration, swing phase duration, stance phase duration, etc.
%
% GENERAL DESCRIPTION
% Enhancing peaks of VT Acc (low-back) by smoothing the signal: using a Gaussian CWT to integrate and differentiate 
% According to: McCamley John, Donati Marco. Grimpampi Eleni, Mazzà Claudia.
% "An enhanced estimate of initial contact and final contact instants of time using lower trunk inertial sensor data."
% Gait and Posture 36,(2012): 316-318"
%
% INPUT
% * Settings. 
%        - Fc_LowPassFilter: cut off frequency used for the low-pass filter
%        - Fc_HighPassFilter: cut-off frequency used for high-pass filters
%        - Fc_HighPassFilter_AutoCorr: cut-off frequency used for high-pass filters, particularly tailored for the autocorrelation signal of the 1st derivative 
%        - MinNumEvents: minimal number of events required to extract temporal gait features
% * AccVT = triaxial acceleration of the defined gait episode (VT axis), downwards (negative) / upwards (positive)
% * Fs = sampling frequency
% * PlotSettings.
%        - StepDetectionDevelopment: 1 (activate plotting), 0 (deactivate)
%       
% HISTORY
% 2019/2ns/September functionized - Encarna Micó Amigo
%
% REMARKS
% * Comments on the right side
% * For peak finding, we considered that the peaks should be more distant for more than half of the dominant frquency


%% [0] DEFINITION OF VARIABLES & PRE-PROCESSING

% 0.1 Definition
Fc_LowPassFilter = Settings.Fc_LowPassFilter;                               % cut off frequency used for the low-pass filter
Fc_HighPassFilter = Settings.Fc_HighPassFilter;                             % cut-off frequency used for high-pass filters
Fc_HighPassFilter_AutoCorr = Settings.Fc_HighPassFilter_AutoCorr;           % cut-off frequency used for high-pass filters, particularly tailored for the autocorrelation signal of the 1st derivative 
GaussianOrder = 10;                                                         % order to calculate wavelet transformation
% 0.2 Pre-Processing
if mean(AccVT) < 0                                                          % signal negative = VT acceleration input is downwards
    AccVT = -AccVT;                                                         % the algorithm requires a positive direction of VT Acc, thus, we ivert the signal (as upwards). If it would be directly provided as upwards, there would not be need to invert it                                                                                        
end
AccVT_Detr = detrend(AccVT);                                                % remove the trend/offset of gravity
AccVT_Detr_Filt = WintFilt_low(AccVT_Detr-mean(AccVT_Detr),Fc_LowPassFilter, Fs); %low-pass filter the signal


%% [1] PEAK ENHANCEMENT: INTEGRATION & DERIVATIVE 

IntegratedSig = cumtrapz(AccVT_Detr_Filt).*1/Fs;                            % integration of signal, which requires to be scaled by 1/Fs
% 1.1 First CWT differentiation  --> initial contact/heel-strike
DerivativeSig1 = derivative_cwt(IntegratedSig','gaus1',GaussianOrder,1./Fs);% transformation with a wavelet function: N vanishing moments ('gaus1' for Gaussian with N=1). % The n-th moment of a function is equal to the n-th derivative of its Fourier transform at zero frequency. So higher the number of zero moments, higher the number of zero derivatives and smoother the signal decays from mid frequency to DC in the frequency domain.
DerivativeSig1_Detrended = detrend(DerivativeSig1);
% 1.2 Second CWT differentiation --> final contact
DerivativeSig2 = derivative_cwt(DerivativeSig1_Detrended,'gaus1',GaussianOrder,1./Fs); % transformation with a 2nd wavelet function
DerivativeSig2_Detrended = detrend(DerivativeSig2);
% 1.3 Plot of obtained signals
if PlotSettings.StepDetectionDevelopment
    figure,plot(AccVT_Detr_Filt,'k'),hold on, plot(DerivativeSig1_Detrended,'r'),plot(DerivativeSig2_Detrended,'b')
end


%% [2.A] PEAK DETECTION - INITIAL CONTACT (IC) / HEEL-STRIKE (HS)           % based on "DerivativeSig1_Detrended"

% 2a.1 Initial peak finding: based on relative amplitude of neighbours
[Amplitude_PeaksFound1,~] = findpeaks(-DerivativeSig1_Detrended); %find the minima of the 1st derivative signal: to avoid detecting "not real peaks"
% 2a.2 Set thresholds - based on amplitude + temporal distance
% - Calculate autocorrelation
% - Extract dominant frequency of autocorrelation signal
% - Take half of the inverse of dominant frequency
ThresholdPeak_Amplitude = mean(Amplitude_PeaksFound1).*0.4;                 % set threshold on a 40% of the average peak amplitude (below that value, in negative = above that value, in positive)
autoCov = xcov(DerivativeSig1_Detrended,'unbiased')';                       % autocorrelation of the signals (evaluation of periodicity)--> function xcov unbiased, for a non attenuate signal.
autoCov_half = autoCov(round(size(autoCov,1)/2)-10:end,1);
autoCov_half_Detr = detrend(autoCov_half);
autoCov_half_Detr_HPFilt = HPFilter(autoCov_half_Detr,Fs,Fc_HighPassFilter_AutoCorr,2);
L = size(autoCov_half_Detr_HPFilt,1);                                       % length of signal
NFFT = 2^nextpow2(L);                                                       % next power of 2 from length of signal
Y = fft(autoCov_half_Detr_HPFilt,NFFT)/L;                                   % Fourier Transform
f = Fs/2*linspace(0,1,NFFT/2+1); warning('off');                            % frequencies of the signal
Power_Spectrum = 2*abs(Y(1:NFFT/2+1));                                      % single-Sided Amplitude Spectrum
[~, Index_Dominantfrequency] = max(Power_Spectrum);                         % max value of Power Spectrum corresponding to the dominant frequency and index
Dominantfrequency = f(Index_Dominantfrequency);
Window = round((1/Dominantfrequency)*Fs);
ThresholdPeak_TempDistance = Window*Settings.PercOfDominantFreq;            % set threshold on half of the dominant frequency (step periodicity) in samples
%ThresholdPeak_TempDistance = 0.3*Fs;                                       % set threshold on a 0.3 seconds of temporal distance (absolute)
% 2a.3 Second peak finding:
[~,InitialContacts] = findpeaks(-DerivativeSig1_Detrended,'MinPeakHeight',ThresholdPeak_Amplitude,'MinPeakDistance',ThresholdPeak_TempDistance);
if PlotSettings.StepDetectionDevelopment
    plot(InitialContacts,DerivativeSig1_Detrended(InitialContacts),'m+')
end


%% [2.B] PEAK DETECTION - FINAL CONTACT (FC)                                % based on "DerivativeSig2_Detrended"

for iIC = 1:size(InitialContacts,2)-1                                       % loop for peak finding between Initial Contact events 
    IC_1st = InitialContacts(iIC);
    IC_2nd = InitialContacts(iIC+1);
    % 2b.1  Initial peak finding: based on relative amplitude of neighbours
    [Amplitude_PeaksFound_BetweenICs1,Sample_PeaksFound_BetweenICs1] = findpeaks(DerivativeSig2_Detrended(IC_1st:IC_2nd)); 
    if PlotSettings.StepDetectionDevelopment
        Junk = DerivativeSig2_Detrended(IC_1st:IC_2nd);
        plot(IC_1st:size(Junk,2)+IC_1st-1,Junk,'g'), hold on,
        plot(Sample_PeaksFound_BetweenICs1+IC_1st-1,Junk(Sample_PeaksFound_BetweenICs1),'c*') % FC
    end
    % 2b.2 Second peak finding: based on a threshold - amplitude
    ThresholdPeak_Amplitude_FC = mean(Amplitude_PeaksFound_BetweenICs1).*0.25; % set threshold on a 25% of the average peak amplitude
    [~,FinalContacts_ShiftedIC] = findpeaks(DerivativeSig2_Detrended(IC_1st:IC_2nd),'MinPeakHeight',ThresholdPeak_Amplitude_FC,'MinPeakDistance',ThresholdPeak_TempDistance);
    if size(FinalContacts_ShiftedIC,2)>1
        [FinalContacts_ShiftedIC,~] = max(FinalContacts_ShiftedIC);         % do not select the first peak, but the highest
    end
    % 2b.3 Add offset (we started search from the first IC)
    if ~isempty(FinalContacts_ShiftedIC)                                    % if no peak with a threshold was found
        FinalContacts(iIC) = IC_1st + FinalContacts_ShiftedIC-1;
    end
end
% 2b.4 Selection of provisional peaks
FinalContacts = nonzeros(FinalContacts)';
if PlotSettings.StepDetectionDevelopment
    plot(FinalContacts,DerivativeSig2_Detrended(FinalContacts),'bs')
end


%% [3] CHECK CORRECT EVENT IDENTIFICATION + CORRECTION

% 3.1 Excluding IC not corresponding to FC
InitialContacts_Selected = [];
FinalContacts_Selected = [];
iStepCount = 1;
for iIC = 1:(min(size(InitialContacts,2),size(FinalContacts,2))-1)          % loop for all events detected, either for all ICs or for all FCs (depending on which has a lower number of events)
    for iFC = 1:size(FinalContacts,2)
        if(FinalContacts(iFC)>InitialContacts(iIC)) && (FinalContacts(iFC)<InitialContacts(iIC+1))
            InitialContacts_Selected(iStepCount)=InitialContacts(iIC);
            FinalContacts_Selected(iStepCount)=FinalContacts(iFC);
            iStepCount = iStepCount+1;
        end
    end
end
% 3.2 Plot detected events
if PlotSettings.StepDetectionDevelopment
    plot(InitialContacts_Selected,DerivativeSig1_Detrended(InitialContacts_Selected),'m*')    
    plot(FinalContacts_Selected,DerivativeSig2_Detrended(FinalContacts_Selected),'bo')        
end


%% [4] CALCULATE DURATION OF CYCLES

IC = InitialContacts_Selected;
FC = FinalContacts_Selected;
if size(IC,2) > Settings.MinNumEvents                                       % at least, for a gait bout counting on xxx steps (number predefined in settings)
    % 4.1 Step durations
    StepDuration_BetweenInitialContacts = diff(IC)'./Fs;                    % Units: seconds ---> scaling by 1/Fs permits to obtain the duration in seconds
    StepDuration_BetweenFinalContacts = diff(FC)'./Fs;                      % Units: seconds
    % 4.2 Stride durations
    for iIC = 1:size(IC,2)-2                                                % -2: because it is for strides, so at least two ICs before the end of the gait (non consecutive, i.e. from the same foot)
        StrideDuration_BetweenInitialContacts(iIC,1) = (IC(iIC+2)-IC(iIC))./Fs; % Units: seconds
        StrideDuration_BetweenFinalContacts(iIC,1) = (FC(iIC+2)-FC(iIC))./Fs; % Units: seconds
    end
    % 4.3 Stance durations
    for iIC = 1:size(IC,2)-1                                                % -1: because it is for stance/swing phases, so at least one IC before the end of the gait
        StanceDuration(iIC,1) = (FC(iIC+1)-IC(iIC))./Fs;                    % Units: seconds
    end
    % 4.4 Swing durations
    SwingDuration = StrideDuration_BetweenInitialContacts - StanceDuration(1:end-1); % Units: seconds
end


%% [5] OUTPUTS

TemporalFeatures.InitialContactSample = IC;
TemporalFeatures.FinalContactsSample = FC;
TemporalFeatures.StepDuration = StepDuration_BetweenInitialContacts;
%TemporalFeatures.StepDuration = StepDuration_BetweenFinalContacts;
TemporalFeatures.StrideDuration = StrideDuration_BetweenInitialContacts;
%TemporalFeatures.StrideDuration = StrideDuration_BetweenFinalContacts;
TemporalFeatures.StanceDuration = StanceDuration;
TemporalFeatures.SwingDuration = SwingDuration;

    
