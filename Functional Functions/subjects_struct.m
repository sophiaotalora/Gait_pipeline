function [data, MetaData] = subjects_struct(Fs, Subjects, DirectoryPath_Data)
    subjectNames = fieldnames(Subjects);
    MetaData = struct;
    MetaData.states=["Test_ON","Test_OFF"];
    data = struct;
    for iSubject = 1%length(subjectNames)
        for iState = 1:length(MetaData.states)
            fullPath = fullfile(DirectoryPath_Data, Subjects.(subjectNames{iSubject}).(MetaData.states{iState}));
            info = h5info(fullPath);
            MetaData.SensorType= {'Right Foot', 'RightFoot'; 'Left Foot', 'LeftFoot'; 'Sternum', 'Sternum'; 'Right Wrist', 'RightWrist'; 'Lumbar', 'LowerBack';'Left Wrist', 'LeftWrist'};
            MetaData.SignalType= {'Acc', 'Accelerometer';  'Gyr', 'Gyroscope'; 'Mag', 'Magnetometer'; 'Timestamp', 'Time'; 'Temp', 'Temperature'};
            MetaData.ImuAcc= {'acc_V','acc_ML','acc_AP'};
            MetaData.ImuGyr={ 'gyro_yaw', 'gyro_pitch','gyro_roll'};                                                          
            MetaData.GroupSize = length(info.Groups(2).Groups);                 %iterate in sensor's groups to create sensor's attributes and codes
            MetaData.SensorSize=length(MetaData.SensorType);
            MetaData.SignalSize= length(MetaData.SignalType);
            MetaData.ImuSize= length(MetaData.ImuAcc);
            
            for iGroup = 1:MetaData.GroupSize
                sensorGroup = info.Groups(2).Groups(iGroup);
                MetaData.attribute(iGroup).value = sensorGroup.Groups(1).Attributes(1).Value;
                MetaData.code(iGroup).value  = strrep(sensorGroup.Name, '/Sensors/', '');
                %MetaData.fs(iGroup).value= sensorGroup.Groups(1).Attributes(12).Value;
            end
            for iSensorType=1:MetaData.SensorSize
                for iSignalType=1:MetaData.SignalSize
                    for icomp=1:MetaData.SensorSize
                        if strcmp(MetaData.SensorType{iSensorType, 1}, MetaData.attribute(icomp).value) 
                            data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).(MetaData.SignalType{iSignalType,1}) = h5read(fullPath, "/Sensors/"+MetaData.code(icomp).value+ "/"+MetaData.SignalType{iSignalType,2})';  
                        end 
                    end %icomp
                                                                       %assign acc and gyr to imu structure (another format for feature extraction)
                    for iImuType=1:MetaData.ImuSize
                        if strcmp(MetaData.SignalType{iSignalType,1}, 'Acc')
                             data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).imu.(MetaData.ImuAcc{iImuType}) =  data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).(MetaData.SignalType{iSignalType,1})(:,iImuType);
                        end
                        if strcmp(MetaData.SignalType{iSignalType,1}, 'Gyr')
                             data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).imu.(MetaData.ImuGyr{iImuType}) =  data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).(MetaData.SignalType{iSignalType,1})(:,iImuType);
                        end
                    end %iImuType
                    if strcmp(MetaData.SignalType{iSignalType}, 'Timestamp')
                        TimeSize= length(data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).(MetaData.SignalType{iSignalType,1}));
                        data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).(MetaData.SignalType{iSignalType,1}) = 0:1/Fs:1/Fs*(TimeSize-1); 
                    end
                end %iSignalType 
                data.("s"+iSubject).(MetaData.states{iState}).SU.(MetaData.SensorType{iSensorType,2}).Fs = 128;
            end %iSensorType                        
        end %iState
    end %iSubject
end %function