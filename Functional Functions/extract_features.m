function data= extract_features(Fs, data, MetaData, Settings, PlotSettings)
    for iSubject = 1%:33%length(subjectNames)
        for iState = 1:length(MetaData.states)           
            WBsize= length(data.("s"+iSubject).(MetaData.states{iState}).SU.WB.straightWalk);
            data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO=struct();  
            Imu =data.("s"+iSubject).(MetaData.states{iState}).SU.LowerBack.imu;
            WB= data.("s"+iSubject).(MetaData.states{iState}).SU.WB.straightWalk;
            for iWB=1:WBsize-1 
                if WB(iWB).start > WB(iWB).end                             % If the start of WB is higuer than the end, then switch
                    Temporal=WB(iWB).start;
                    WB(iWB).start = WB(iWB).end;
                    WB(iWB).end = Temporal;
                end
                if (WB(iWB).end - WB(iWB).start) < Fs*2                    % If it is less than 2 seconds the WB, then skip to next one
                    iWB= iWB+1;
                    if iWB >= WBsize
                        continue;
                    end
                end

                % Pace
                Stepres= extract_StepLengthVelocityMean_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepLength_Mean(iWB) = Stepres.StepLength_Mean;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepVelocity_Mean(iWB) = Stepres.StepVelocity_Mean;

                %% Rythm
                %data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.cadence(iWB) = extract_Cadence_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                Strideres = extract_StepStrideDurationMean_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepDuration_Mean(iWB) = Strideres.StepDuration_Mean;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideDuration_Mean(iWB) = Strideres.StrideDuration_Mean;
                
                Swingstance = extract_SwingStanceDurationMean_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDuration_Mean(iWB)= Swingstance.SwingDuration_Mean;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDurationPercentage_Mean(iWB)= Swingstance.SwingDurationPercentage_Mean;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDuration_Mean(iWB)= Swingstance.StanceDuration_Mean;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDurationPercentage_Mean(iWB)= Swingstance.StanceDurationPercentage_Mean;
                
                %% Variability
                Steplen = extract_StepLengthVelocityVariability_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepLength_CoeffVariation(iWB) = Steplen.StepLength_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepLength_SD(iWB) = Steplen.StepLength_SD;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepVelocity_CoeffVariation(iWB) = Steplen.StepVelocity_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepVelocity_SD(iWB) = Steplen.StepVelocity_SD;
                Stepdurres= extract_StepStrideDurationVariability_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepDuration_CoeffVariation(iWB) = Stepdurres.StepDuration_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepDuration_SD(iWB) = Stepdurres.StepDuration_SD;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideDuration_CoeffVariation(iWB) = Stepdurres.StrideDuration_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideDuration_SD(iWB) = Stepdurres.StrideDuration_SD;
                
                % Swing
                Swingdurres= extract_SwingStanceDurationVariability_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); 
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDuration_CoeffVariation(iWB) = Swingdurres.SwingDuration_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDurationPercentage_CoeffVariation(iWB) = Swingdurres.SwingDurationPercentage_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDuration_SD(iWB) = Swingdurres.SwingDuration_SD;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDurationPercentage_SD(iWB) = Swingdurres.SwingDurationPercentage_SD;
                % Stance
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDuration_CoeffVariation(iWB) = Swingdurres.StanceDuration_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDurationPercentage_CoeffVariation(iWB) = Swingdurres.StanceDurationPercentage_CoeffVariation;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDuration_SD(iWB) = Swingdurres.StanceDuration_SD;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDurationPercentage_SD(iWB) = Swingdurres.StanceDurationPercentage_SD;
                
                
                %% Symmetry
                Stepassym = extract_StepLengthVelocityAsymmetry_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepLength_Asymmetry(iWB) = Stepassym.StepLength_Asymmetry;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepVelocity_Asymmetry(iWB) = Stepassym.StepVelocity_Asymmetry;
                Stepdurassym = extract_StepStrideDurationAsymmetry_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepDuration_Asymmetry(iWB) = Stepdurassym.StepDuration_Asymmetry;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideDuration_Asymmetry(iWB) = Stepdurassym.StrideDuration_Asymmetry;
                Swingdurassym = extract_SwingStanceDurationAsymmetry_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDuration_Asymmetry(iWB) = Swingdurassym.SwingDuration_Asymmetry;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDuration_Asymmetry(iWB) = Swingdurassym.StanceDuration_Asymmetry;
                
                %% Regularity
                Regularityacc= extract_Regularity_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepRegularity_Acc_VT(iWB)=Regularityacc.StepRegularity_Acc_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepRegularity_Acc_AP(iWB)=Regularityacc.StepRegularity_Acc_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideRegularity_Acc_VT(iWB)=Regularityacc.StrideRegularity_Acc_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideRegularity_Acc_AP(iWB)=Regularityacc.StrideRegularity_Acc_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationRatio_VT(iWB)=Regularityacc.Symmetry_AutocorrelationRatio_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationRatio_AP(iWB)=Regularityacc.Symmetry_AutocorrelationRatio_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationDiff_VT(iWB)=Regularityacc.Symmetry_AutocorrelationDiff_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationDiff_AP(iWB)=Regularityacc.Symmetry_AutocorrelationDiff_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.GaitSymmetryIndex_Autocorr(iWB) = extract_Regularity_GaitSymmetryIndex_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                
                %% Spectral Components
                Spectral_Acc = extract_contentSpectralPowerDensity_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                Spectral_Gyr = extract_contentSpectralPowerDensity_Gyro(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                HarmonicRatio_Acc = extract_HarmonicRatio_PerStrides_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                HarmonicRatio_Gyr = extract_HarmonicRatio_PerStrides_Gyro(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                NameChannelsAcc = {'VT','ML','AP'};
                NameChannelsGyr = {'Yaw','Pitch','Roll'};
                VariablesPhaseNames = {'DominantFreq','Amplitude','AmplitudeNorm','Width','Slope','Range','MeanPower','MedianPower','MeanFreq','MedianFreq','WidthNorm','SlopeNorm','RangeNorm','IntegratedPower','HarmonicRatio','IndexHarmonicity'};
                for iChannel = 1:size(NameChannelsAcc,2)
                    data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.("HarmonicRatioAcc_"+(NameChannelsAcc{iChannel}))(iWB)=HarmonicRatio_Acc.HarmonicRatio_AvPerStrides_Acc.(NameChannelsAcc{iChannel});
                    for iFeature = 1:size(VariablesPhaseNames,2)
                        data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.("SpectralAcc_"+(NameChannelsAcc{iChannel})+"_"+(VariablesPhaseNames{iFeature}))(iWB)=Spectral_Acc.SpectralFeatures_Acc.(NameChannelsAcc{iChannel}).(VariablesPhaseNames{iFeature});  
                    end
                end
                for iChannel = 1:size(NameChannelsGyr,2)
                    data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.("HarmonicRatioGyr_"+(NameChannelsGyr{iChannel}))(iWB)=HarmonicRatio_Gyr.HarmonicRatio_AvPerStrides_Gyr.(NameChannelsGyr{iChannel});
                    for iFeature = 1:size(VariablesPhaseNames,2)
                        data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.("SpectralGyr_"+(NameChannelsGyr{iChannel})+"_"+(VariablesPhaseNames{iFeature}))(iWB)=Spectral_Gyr.SpectralFeatures_Gyr.(NameChannelsGyr{iChannel}).(VariablesPhaseNames{iFeature});
                    end
                end
                
                %% MAGNITUDE
                % % Acceleration
                RMS_Acc= extract_RMS_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);   % RMS pure of acceleration
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Acc_VT(iWB)=RMS_Acc.RMSq_SigComplete.VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Acc_ML(iWB)=RMS_Acc.RMSq_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Acc_AP(iWB)=RMS_Acc.RMSq_SigComplete.AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Acc_Combined(iWB)=RMS_Acc.RMSq_SigComplete.Combined;

                RMSRatio_Acc = extract_RMSRatio_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % RMS ratio of acceleration (RMS of each component / RMS of resultant);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMSRatio_Acc_VT(iWB)=RMSRatio_Acc.RMSratio_SigComplete.VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMSRatio_Acc_ML(iWB)=RMSRatio_Acc.RMSratio_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMSRatio_Acc_AP(iWB)=RMSRatio_Acc.RMSratio_SigComplete.AP;
               
                % % Angular Velocity
          
                RMS_Gyro = extract_RMS_Gyro(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);  % RMS pure of angular velocity
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Gyro_Yaw(iWB)=RMS_Gyro.RMSq_SigComplete.Yaw;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Gyro_Pitch(iWB)=RMS_Gyro.RMSq_SigComplete.Pitch;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Gyro_Roll(iWB)=RMS_Gyro.RMSq_SigComplete.Roll;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMS_Gyro_Combined(iWB)=RMS_Gyro.RMSq_SigComplete.Combined;
           
                RMSRatio_Gyro = extract_RMSRatio_Gyro(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % RMS ratio of angular velocity (RMS of each component / RMS of resultant)
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMSRatio_Gyro_Yaw(iWB)=RMSRatio_Gyro.RMSratio_SigComplete.Yaw;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMSRatio_Gyro_Pitch(iWB)=RMSRatio_Gyro.RMSratio_SigComplete.Pitch;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.RMSRatio_Gyro_Roll(iWB)=RMSRatio_Gyro.RMSratio_SigComplete.Roll;
                
                %% COMPLEXITY
                % % Jerk: 1st derivative of acceleration
               
                JerkRMS_Acc = extract_JerkRMS_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % RMS of jerk
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRMS_Acc_VT(iWB)=JerkRMS_Acc.JerkRMS_SigComplete.VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRMS_Acc_ML(iWB)=JerkRMS_Acc.JerkRMS_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRMS_Acc_AP(iWB)=JerkRMS_Acc.JerkRMS_SigComplete.AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRMS_Acc_Combined(iWB)=JerkRMS_Acc.JerkRMS_SigComplete.Combined;
  
                JerkMeanLogRatio_Acc = extract_JerkMeanLogRatio_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % Ratio of jerk, i.e. jerk on each component / jerk on resultant
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMeanLogRatio_Acc_ML(iWB)=JerkMeanLogRatio_Acc.JerkMeanLogRatio_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMeanLogRatio_Acc_AP(iWB)=JerkMeanLogRatio_Acc.JerkMeanLogRatio_SigComplete.AP;
                         
                JerkMax_Acc = extract_JerkMax_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % RMS of jerk
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMax_Acc_VT(iWB)=JerkMax_Acc.JerkMax_SigComplete.VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMax_Acc_ML(iWB)=JerkMax_Acc.JerkMax_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMax_Acc_AP(iWB)=JerkMax_Acc.JerkMax_SigComplete.AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMax_Acc_Combined(iWB)=JerkMax_Acc.JerkMax_SigComplete.Combined;
            
                JerkMin_Acc = extract_JerkMin_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % RMS of jerk
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMin_Acc_VT(iWB)=JerkMin_Acc.JerkMin_SigComplete.VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMin_Acc_ML(iWB)=JerkMin_Acc.JerkMin_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMin_Acc_AP(iWB)=JerkMin_Acc.JerkMin_SigComplete.AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkMin_Acc_Combined(iWB)=JerkMin_Acc.JerkMin_SigComplete.Combined;
           
                JerkRange_Acc = extract_JerkRange_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % Min of jerk
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRange_Acc_VT(iWB)=JerkRange_Acc.JerkRange_SigComplete.VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRange_Acc_ML(iWB)=JerkRange_Acc.JerkRange_SigComplete.ML;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRange_Acc_AP(iWB)=JerkRange_Acc.JerkRange_SigComplete.AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.JerkRange_Acc_Combined(iWB)=JerkRange_Acc.JerkRange_SigComplete.Combined;
 
                % % Coefficient of attenuation: we need data from sensors located at two different locations
                % % Lyapunov exponent
                Types=["LyapunovRC_Acc", "LyapunovW_Acc"];
              
                LyapunovExponent_Acc = extract_LyapunovExponent_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); 
               
                for iType= 1:length(Types)
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_VT")(iWB)=LyapunovExponent_Acc.(Types(iType)).VT;
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_ML")(iWB)=LyapunovExponent_Acc.(Types(iType)).ML;
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_AP")(iWB)=LyapunovExponent_Acc.(Types(iType)).AP;
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_Combined")(iWB)=LyapunovExponent_Acc.(Types(iType)).Combined;
                end
                
                %Gyro
                Types=["LyapunovRC_Gyr", "LyapunovW_Gyr"];
                
                LyapunovExponent_Gyro = extract_LyapunovExponent_Gyro(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); 
               
                for iType= 1:length(Types)
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_Yaw")(iWB)=LyapunovExponent_Gyro.(Types(iType)).Yaw;
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_Pitch")(iWB)=LyapunovExponent_Gyro.(Types(iType)).Pitch;
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_Roll")(iWB)=LyapunovExponent_Gyro.(Types(iType)).Roll;
                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.((Types(iType))+"_Combined")(iWB)=LyapunovExponent_Gyro.(Types(iType)).Combined;
                end
                
                % % Correlation between accelerometry and angular velocity
               
                Correlation_AccGyro = extract_correlation_AccGyro(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); 
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Correlation_AccGyroVT_Yaw(iWB)=Correlation_AccGyro.CorrAccGyr_SigComplete.VT_Yaw;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Correlation_AccGyroML_Pitch(iWB)=Correlation_AccGyro.CorrAccGyr_SigComplete.ML_Pitch;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Correlation_AccGyro_AP_Roll(iWB)=Correlation_AccGyro.CorrAccGyr_SigComplete.AP_Roll;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Correlation_AccGyro_Combined2(iWB)=Correlation_AccGyro.CorrAccGyr_SigComplete.Combined2;
               
                %}
                
            end %iWB       
        end %iState
    end %iSubject
    toc
end