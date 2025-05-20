function data= extract_features_step(Fs, data, Subjects, MetaData, Settings, PlotSettings)
    for iSubject = 1:33%:length(subjectNames)
        for iState = 1:length(MetaData.states)           
            WBsize= length(data.("s"+iSubject).(MetaData.states{iState}).SU.WB.straightWalk);
            StepSize= length(data.("s"+iSubject).(MetaData.states{iState}).SU.WB.FirstSteps);
            data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO=struct();  
            Imu =data.("s"+iSubject).(MetaData.states{iState}).SU.LowerBack.imu;
            WB= data.("s"+iSubject).(MetaData.states{iState}).SU.WB.straightWalk;
            Steps= data.("s"+iSubject).(MetaData.states{iState}).SU.WB.FirstSteps;
            for iWB= 1:WBsize-1 
                Sections= ["First", "Middle", "Last"];
                %{
                % Pace
                try
                    Stepres= extract_StepLengthVelocityMean_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                    Names = ["StepLength","StepVelocity"];
                    for iName = 1:length(Names)
                        for iSec = 1:length(Sections)
                            SecName= char(Sections(iSec));
                            VarName = Names{iName};
                            data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.(SecName).(VarName)(:,iWB) = Stepres.(SecName).(VarName);
                        end
                    end
%                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.First.StepLength_Mean(:,iWB) = Stepres.First.StepLength;
%                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Middle.StepLength_Mean(:,iWB) = Stepres.Middle.StepLength;
%                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Last.StepLength_Mean(:,iWB) = Stepres.Last.StepLength;
%                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.First.StepVelocity_Mean(:,iWB) = Stepres.First.StepVelocity;
%                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Middle.StepVelocity_Mean(:,iWB) = Stepres.Middle.StepVelocity;
%                     data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Last.StepVelocity_Mean(:,iWB) = Stepres.Last.StepVelocity;
                catch ME
                    disp('An error occurred, skipping to the next WB');
                    %disp(ME.message)
                    %disp(iWB)
                    %disp(iSubject)
                    iWB = iWB + 1;
                    if iWB >= length(iWB)
                        continue;
                    else
                    Stepres= extract_StepLengthVelocityMean_Acc(Imu,Fs,WB(iWB,:),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                    Names = ["StepLength","StepVelocity"];
                    for iName = 1:length(Names)
                        for iSec = 1:length(Sections)
                            SecName= char(Sections(iSec));
                            VarName = Names{iName};
                            data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.(SecName).(VarName)(:,iWB) = Stepres.(SecName).(VarName);
                        end
                    end
                    end
                end
                
                % Rythm
                %data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.cadence(iWB) = extract_Cadence_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                
                Strideres = extract_StepStrideDurationMean_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                Names =  ["StepDuration","StrideDuration"];
                for iName = 1:length(Names)
                    for iSec = 1:length(Sections)
                        SecName= char(Sections(iSec));
                        VarName = Names{iName};
                        data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.(SecName).(VarName)(:,iWB) = Strideres.(SecName).(VarName);
                    end
                end
               
                
                Swingstance = extract_SwingStanceDurationMean_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                Names= ["SwingDuration","SwingDurationPercentage", "StanceDuration","StanceDurationPercentage"];
                for iName = 1:length(Names)
                    for iSec = 1:length(Sections)
                        SecName= char(Sections(iSec));
                        VarName = Names{iName};
                        data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.(SecName).(VarName)(:,iWB) = Swingstance.(SecName).(VarName);
                    end
                end  
               
               
                % Variability
                Steplen = extract_StepLengthVelocityVariability_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                
                Names= ["StepLength_CoeffVariation","StepLength_SD", "StepVelocity_CoeffVariation","StepVelocity_SD"];
                for iName = 1:length(Names)
                    for iSec = 1:length(Sections)
                        SecName= char(Sections(iSec));
                        VarName = Names{iName};
                        data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.(SecName).(VarName)(:,iWB) = Steplen.(SecName).(VarName);
                    end
                end   
                %}
                try
                    Stepdurres= extract_StepStrideDurationVariability_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                    Names= ["StepDuration_CoeffVariation","StepDuration_SD", "StrideDuration_CoeffVariation","StrideDuration_SD"];
                    for iName = 1:length(Names)
                        for iSec = 1:length(Sections)
                            SecName= char(Sections(iSec));
                            VarName = Names{iName};
                            data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.(SecName).(VarName)(:,iWB) = Stepdurres.(SecName).(VarName);
                        end
                    end  
                catch ME
                    iWB = iWB + 1;
                    if iWB >= length(iWB)
                        continue;
                    end
                end
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
                %{
                % Symmetry
                Stepassym = extract_StepLengthVelocityAsymmetry_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepLength_Asymmetry(iWB) = Stepassym.StepLength_Asymmetry;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepVelocity_Asymmetry(iWB) = Stepassym.StepVelocity_Asymmetry;
                Stepdurassym = extract_StepStrideDurationAsymmetry_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepDuration_Asymmetry(iWB) = Stepdurassym.StepDuration_Asymmetry;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideDuration_Asymmetry(iWB) = Stepdurassym.StrideDuration_Asymmetry;
                Swingdurassym = extract_SwingStanceDurationAsymmetry_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.SwingDuration_Asymmetry(iWB) = Swingdurassym.SwingDuration_Asymmetry;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StanceDuration_Asymmetry(iWB) = Swingdurassym.StanceDuration_Asymmetry;
                
                %Regularity
                Regularityacc= extract_Regularity_Acc(Imu,Fs,Steps(iSteps),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepRegularity_Acc_VT(iWB)=Regularityacc.StepRegularity_Acc_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StepRegularity_Acc_AP(iWB)=Regularityacc.StepRegularity_Acc_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideRegularity_Acc_VT(iWB)=Regularityacc.StrideRegularity_Acc_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.StrideRegularity_Acc_AP(iWB)=Regularityacc.StrideRegularity_Acc_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationRatio_VT(iWB)=Regularityacc.Symmetry_AutocorrelationRatio_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationRatio_AP(iWB)=Regularityacc.Symmetry_AutocorrelationRatio_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationDiff_VT(iWB)=Regularityacc.Symmetry_AutocorrelationDiff_VT;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.Symmetry_AutocorrelationDiff_AP(iWB)=Regularityacc.Symmetry_AutocorrelationDiff_AP;
                data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO.GaitSymmetryIndex_Autocorr(iWB) = extract_Regularity_GaitSymmetryIndex_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings);
                %}
                %{
                %Spectral Components
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
                %}
                %Complexity
                %JerkRMS_Acc = extract_JerkRMS_Acc(Imu,Fs,WB(iWB),data.("s"+iSubject).(MetaData.states{iState}).SU.WB_SDMO,Settings,PlotSettings); % RMS of jerk
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
        end %iState
    end %iSubject
    toc
end