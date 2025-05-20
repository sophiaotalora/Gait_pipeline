function subjects = processfilesh5(DirectoryPath_Data)
    folderPath = DirectoryPath_Data;
    files = dir(fullfile(folderPath, '*.h5'));
    subjects = struct();
    subjectCounter = 1;

    for i = 1:length(files)
        fileName = files(i).name;
        dateStr = fileName(1:8);                                           %date and time from the file name
        timeStr = fileName(10:15); 
                                                                           %datetime format
        fileDateTime = datetime([dateStr, timeStr], 'InputFormat', 'yyyyMMddHHmmss');
        subjectID = findSubjectID(subjects, dateStr);                      %date exists in the structure?
        if isempty(subjectID)                                              %create new subject ID
            subjectID = ['s', num2str(subjectCounter)];
            subjectCounter = subjectCounter + 1;
        end

        if ~isfield(subjects, subjectID)                                   %subject's file list
            subjects.(subjectID).files = struct('fileName', fileName, 'dateTime', fileDateTime);
        else
            subjects.(subjectID).files(end+1) = struct('fileName', fileName, 'dateTime', fileDateTime);
        end
    end

    %"ON" and "OFF" status for each subject based on time
    subjectNames = fieldnames(subjects);
    for i = 1:length(subjectNames)
        subjectFiles = subjects.(subjectNames{i}).files;
        allTimes = [subjectFiles.dateTime];
        
        % Find the earliest and latest times
        [earliestTime, earliestIdx] = min(allTimes);
        [latestTime, latestIdx] = max(allTimes);
        
        % Assign "ON" and "OFF" based on times
        if date(earliestTime) == 
            % If the times are on the same date
        subjects.(subjectNames{i}).Test_ON = subjectFiles(latestIdx).fileName;
        subjects.(subjectNames{i}).Test_OFF = subjectFiles(earliestIdx).fileName;
        
         

        % Clean up the temporary 'files' field
        subjects.(subjectNames{i}) = rmfield(subjects.(subjectNames{i}), 'files');
    end
end

%find the subject ID based on the date
function subjectID = findSubjectID(subjects, dateStr)
    subjectID = '';  
    subjectNames = fieldnames(subjects);
    for i = 1:length(subjectNames)
        if strcmp(subjects.(subjectNames{i}).files(1).fileName(1:8), dateStr)
            subjectID = subjectNames{i};
            return;
        end
    end
end