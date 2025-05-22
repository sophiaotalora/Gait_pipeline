function Data= extract_features_c(Fs, Data, Settings, PlotSettings)
    States = fieldnamesr(Data, 2);
    States = unique(regexprep(States, '.*\.', ''));
    SubjectNames = fieldnames(Data);
    Trials= fieldnamesr(Data, 3);
    Trials = unique(regexprep(Trials, '.*\.', ''));
    for iSubject = 1:length(SubjectNames)
        for iState = 1:length(States)
            for iTrial = 1:length(Trials)         
                WBsize= length(Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.straightWalk);
                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO=struct();  
                Imu =Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).imu;
                WB= Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.straightWalk;
                for iWB=1:WBsize-1                    
                    % Pace
                    try
                        Stepres= extract_StepLengthVelocityMean_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepLength_Mean(iWB) = Stepres.StepLength_Mean;
                        Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepVelocity_Mean(iWB) = Stepres.StepVelocity_Mean;
                    catch ME
                        disp(ME.message);
                        %disp('An error occurred, skipping to the next iWB');
                        iWB = iWB + 1;
                        if iWB >= length(iWB)
                            continue;
                        else
                        Stepres= extract_StepLengthVelocityMean_Acc(Imu,Fs,WB(iWB),Data.Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepLength_Mean(iWB) = Stepres.StepLength_Mean;
                        Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepVelocity_Mean(iWB) = Stepres.StepVelocity_Mean;
                        end
                    end
                    %
                    % Rythm
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.cadence(iWB) = extract_Cadence_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Strideres = extract_StepStrideDurationMean_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepDuration_Mean(iWB) = Strideres.StepDuration_Mean;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StrideDuration_Mean(iWB) = Strideres.StrideDuration_Mean;
                    
                    Swingstance = extract_SwingStanceDurationMean_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDuration_Mean(iWB)= Swingstance.SwingDuration_Mean;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDurationPercentage_Mean(iWB)= Swingstance.SwingDurationPercentage_Mean;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDuration_Mean(iWB)= Swingstance.StanceDuration_Mean;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDurationPercentage_Mean(iWB)= Swingstance.StanceDurationPercentage_Mean;
                    
                    % Variability
                    Steplen = extract_StepLengthVelocityVariability_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepLength_CoeffVariation(iWB) = Steplen.StepLength_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepLength_SD(iWB) = Steplen.StepLength_SD;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepVelocity_CoeffVariation(iWB) = Steplen.StepVelocity_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepVelocity_SD(iWB) = Steplen.StepVelocity_SD;
                    Stepdurres= extract_StepStrideDurationVariability_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepDuration_CoeffVariation(iWB) = Stepdurres.StepDuration_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepDuration_SD(iWB) = Stepdurres.StepDuration_SD;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StrideDuration_CoeffVariation(iWB) = Stepdurres.StrideDuration_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StrideDuration_SD(iWB) = Stepdurres.StrideDuration_SD;
                    % Swing
                    Swingdurres= extract_SwingStanceDurationVariability_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings); 
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDuration_CoeffVariation(iWB) = Swingdurres.SwingDuration_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDurationPercentage_CoeffVariation(iWB) = Swingdurres.SwingDurationPercentage_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDuration_SD(iWB) = Swingdurres.SwingDuration_SD;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDurationPercentage_SD(iWB) = Swingdurres.SwingDurationPercentage_SD;
                    % Stance
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDuration_CoeffVariation(iWB) = Swingdurres.StanceDuration_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDurationPercentage_CoeffVariation(iWB) = Swingdurres.StanceDurationPercentage_CoeffVariation;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDuration_SD(iWB) = Swingdurres.StanceDuration_SD;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDurationPercentage_SD(iWB) = Swingdurres.StanceDurationPercentage_SD;
                    
                    % Symmetry
                    Stepassym = extract_StepLengthVelocityAsymmetry_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepLength_Asymmetry(iWB) = Stepassym.StepLength_Asymmetry;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepVelocity_Asymmetry(iWB) = Stepassym.StepVelocity_Asymmetry;
                    Stepdurassym = extract_StepStrideDurationAsymmetry_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepDuration_Asymmetry(iWB) = Stepdurassym.StepDuration_Asymmetry;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StrideDuration_Asymmetry(iWB) = Stepdurassym.StrideDuration_Asymmetry;
                    Swingdurassym = extract_SwingStanceDurationAsymmetry_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.SwingDuration_Asymmetry(iWB) = Swingdurassym.SwingDuration_Asymmetry;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StanceDuration_Asymmetry(iWB) = Swingdurassym.StanceDuration_Asymmetry;
                    %
                    %Regularity
                    Regularityacc= extract_Regularity_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepRegularity_Acc_VT(iWB)=Regularityacc.StepRegularity_Acc_VT;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StepRegularity_Acc_AP(iWB)=Regularityacc.StepRegularity_Acc_AP;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StrideRegularity_Acc_VT(iWB)=Regularityacc.StrideRegularity_Acc_VT;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.StrideRegularity_Acc_AP(iWB)=Regularityacc.StrideRegularity_Acc_AP;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.Symmetry_AutocorrelationRatio_VT(iWB)=Regularityacc.Symmetry_AutocorrelationRatio_VT;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.Symmetry_AutocorrelationRatio_AP(iWB)=Regularityacc.Symmetry_AutocorrelationRatio_AP;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.Symmetry_AutocorrelationDiff_VT(iWB)=Regularityacc.Symmetry_AutocorrelationDiff_VT;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.Symmetry_AutocorrelationDiff_AP(iWB)=Regularityacc.Symmetry_AutocorrelationDiff_AP;
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.GaitSymmetryIndex_Autocorr(iWB) = extract_Regularity_GaitSymmetryIndex_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                    try
                        %Spectral Components
                        Spectral_Acc = extract_contentSpectralPowerDensity_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        Spectral_Gyr = extract_contentSpectralPowerDensity_Gyro(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        HarmonicRatio_Acc = extract_HarmonicRatio_PerStrides_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        HarmonicRatio_Gyr = extract_HarmonicRatio_PerStrides_Gyro(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        NameChannelsAcc = {'VT','ML','AP'};
                        NameChannelsGyr = {'Yaw','Pitch','Roll'};
                        VariablesPhaseNames = {'DominantFreq','Amplitude','AmplitudeNorm','Width','Slope','Range','MeanPower','MedianPower','MeanFreq','MedianFreq','WidthNorm','SlopeNorm','RangeNorm','IntegratedPower','HarmonicRatio','IndexHarmonicity'};
                        for iChannel = 1:size(NameChannelsAcc,2)
                            Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("HarmonicRatioAcc_"+(NameChannelsAcc{iChannel}))(iWB)=HarmonicRatio_Acc.HarmonicRatio_AvPerStrides_Acc.(NameChannelsAcc{iChannel});
                            for iFeature = 1:size(VariablesPhaseNames,2)
                                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("SpectralAcc_"+(NameChannelsAcc{iChannel})+"_"+(VariablesPhaseNames{iFeature}))(iWB)=Spectral_Acc.SpectralFeatures_Acc.(NameChannelsAcc{iChannel}).(VariablesPhaseNames{iFeature});  
                            end
                        end
                        for iChannel = 1:size(NameChannelsGyr,2)
                            Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("HarmonicRatioGyr_"+(NameChannelsGyr{iChannel}))(iWB)=HarmonicRatio_Gyr.HarmonicRatio_AvPerStrides_Gyr.(NameChannelsGyr{iChannel});
                            for iFeature = 1:size(VariablesPhaseNames,2)
                                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("SpectralGyr_"+(NameChannelsGyr{iChannel})+"_"+(VariablesPhaseNames{iFeature}))(iWB)=Spectral_Gyr.SpectralFeatures_Gyr.(NameChannelsGyr{iChannel}).(VariablesPhaseNames{iFeature});
                            end
                        end
                    catch 
                        iWB = iWB + 1;
                        if iWB >= length(iWB)
                            continue;
                        else
                            Spectral_Acc = extract_contentSpectralPowerDensity_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        Spectral_Gyr = extract_contentSpectralPowerDensity_Gyro(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        HarmonicRatio_Acc = extract_HarmonicRatio_PerStrides_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        HarmonicRatio_Gyr = extract_HarmonicRatio_PerStrides_Gyro(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings);
                        NameChannelsAcc = {'VT','ML','AP'};
                        NameChannelsGyr = {'Yaw','Pitch','Roll'};
                        VariablesPhaseNames = {'DominantFreq','Amplitude','AmplitudeNorm','Width','Slope','Range','MeanPower','MedianPower','MeanFreq','MedianFreq','WidthNorm','SlopeNorm','RangeNorm','IntegratedPower','HarmonicRatio','IndexHarmonicity'};
                        for iChannel = 1:size(NameChannelsAcc,2)
                            Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("HarmonicRatioAcc_"+(NameChannelsAcc{iChannel}))(iWB)=HarmonicRatio_Acc.HarmonicRatio_AvPerStrides_Acc.(NameChannelsAcc{iChannel});
                            for iFeature = 1:size(VariablesPhaseNames,2)
                                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("SpectralAcc_"+(NameChannelsAcc{iChannel})+"_"+(VariablesPhaseNames{iFeature}))(iWB)=Spectral_Acc.SpectralFeatures_Acc.(NameChannelsAcc{iChannel}).(VariablesPhaseNames{iFeature});  
                            end
                        end
                        for iChannel = 1:size(NameChannelsGyr,2)
                            Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("HarmonicRatioGyr_"+(NameChannelsGyr{iChannel}))(iWB)=HarmonicRatio_Gyr.HarmonicRatio_AvPerStrides_Gyr.(NameChannelsGyr{iChannel});
                            for iFeature = 1:size(VariablesPhaseNames,2)
                                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO.("SpectralGyr_"+(NameChannelsGyr{iChannel})+"_"+(VariablesPhaseNames{iFeature}))(iWB)=Spectral_Gyr.SpectralFeatures_Gyr.(NameChannelsGyr{iChannel}).(VariablesPhaseNames{iFeature});
                            end
                        end
                        end
                    end
                    %}
                    %Complexity
                    %JerkRMS_Acc = extract_JerkRMS_Acc(Imu,Fs,WB(iWB),Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB_SDMO,Settings,PlotSettings); % RMS of jerk
    
                    %}
         
                    %% MAGNITUDE
                    %% Acceleration
                    % [WB_SDMO] = extract_RMS_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings);   % RMS pure of acceleration
                    % [WB_SDMO] = extract_RMSRatio_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % RMS ratio of acceleration (RMS of each component / RMS of resultant)
                    % % Angular Velocity
                    % [WB_SDMO] = extract_RMS_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings);  % RMS pure of angular velocity
                    % [WB_SDMO] = extract_RMSRatio_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % RMS ratio of angular velocity (RMS of each component / RMS of resultant)
                    % 
                    %  m 
                    %% COMPLEXITY
                    % % Jerk: 1st derivative of acceleration
                    % [WB_SDMO] = extract_JerkRMS_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % RMS of jerk
                    % [WB_SDMO] = extract_JerkMeanLogRatio_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % Ratio of jerk, i.e. jerk on each component / jerk on resultant
                    % [WB_SDMO] = extract_JerkMax_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % Max of jerk
                    % [WB_SDMO] = extract_JerkMin_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % Min of jerk
                    % [WB_SDMO] = extract_JerkRange_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings); % Range of jerk
                    % % Coefficient of attenuation: we need data from sensors located at two different locations
                    %% Lyapunov exponent
                    % [WB_SDMO] = extract_LyapunovExponent_Acc(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings);
                    % [WB_SDMO] = extract_LyapunovExponent_Gyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings);
                    % % Correlation between accelerometry and angular velocity
                    % [WB_SDMO] = extract_correlation_AccGyro(imu,fs,WB_WBD,WB_SDMO,Settings,PlotSettings);                
                end %iWB
            end %iTrial
        end %iState
    end %iSubject