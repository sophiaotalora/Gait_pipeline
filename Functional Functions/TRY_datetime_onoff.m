clc;clear all; close all;
folderPath = 'C:\Users\sophi\Heriot-Watt University\Micó Amigo, Encarna - Digital Mobility Technology\Data\OHSU - PD medication\Standardized Data\Raw_data';

files = dir(fullfile(folderPath, '*.h5'));                                 %list of all the .h5 files in the folder
subjects = struct();
subjectCounter = 1;
for i = 1:length(files)    
    fileName = files(i).name;    
    dateStr = fileName(1:8);                                               %date and time from the file name
    timeStr = fileName(10:15);  
    fileDateTime = datetime([dateStr, timeStr], 'InputFormat', 'yyyyMMddHHmmss'); %datetime format

    subjectID = findSubjectID(subjects, dateStr);                          %verify if date already exists in the structure   
    if isempty(subjectID)
        subjectID = ['s', num2str(subjectCounter)];                        %create new subject ID
        subjectCounter = subjectCounter + 1;
    end
    filePath = fullfile(folderPath, fileName);
    fileData = h5read(filePath, '/');                                      %load all data from the .h5 file

    if ~isfield(subjects, subjectID)                                       %append file data to subject's file list
        subjects.(subjectID).files = struct('dateTime', fileDateTime, 'date', fileData);
    else
        subjects.(subjectID).files(end+1) = struct('dateTime', fileDateTime, 'data', fileData);
    end
end
%%
%"ON" and "OFF" status for each subject based on time
subjectNames = fieldnames(subjects);
%length(subjectNames)
%hacer load directamente en el file - validation cadence (bloque 2 - loading imus, open y load)
for i = 1:length(subjectNames)
    subjectFiles = subjects.(subjectNames{i}).files;                       %list of files and their datetimes
    allTimes = [subjectFiles.dateTime];    
    [~, earliestIdx] = min(allTimes);                                      %find the earliest and latest times
    [~, latestIdx] = max(allTimes);
    cd(folderPath);
    subjects.(subjectNames{i}).ON = subjectFiles(latestIdx).fileName;      %"ON" to the latest time and "OFF" to the earliest time
    subjects.(subjectNames{i}).OFF = subjectFiles(earliestIdx).fileName;
    
    subjects.(subjectNames{i}) = rmfield(subjects.(subjectNames{i}), 'files');  % Remove 'files' field
end


%find the subject ID based on the date
function subjectID = findSubjectID(subjects, dateStr)                                                                        
    subjectID = '';                                                        %subjectID as empty                                                    
    subjectNames = fieldnames(subjects);
    for i = 1:length(subjectNames)                                         %Iterate over the existing fields in subjects
        if strcmp(subjects.(subjectNames{i}).files(1).fileName(1:8), dateStr) %Compare the dates
            subjectID = subjectNames{i};
            return;
        end
    end
end
%%