function subjects = processFiles(DirectoryPath_Data)
    % Directory containing the .h5 files
    folderPath = DirectoryPath_Data;

    % Get a list of all the .h5 files in the folder
    files = dir(fullfile(folderPath, '*.h5'));
    
    % Initialize an empty structure array to store subject data
    subjects = struct();
    subjectCounter = 1;

    % Iterate over each file
    for i = 1:length(files)
        % Get the file name
        fileName = files(i).name;

        % Extract the date and time from the file name
        dateStr = fileName(1:8);  % 'yyyyMMdd'
        timeStr = fileName(10:15); % 'HHmmss'

        % Convert date and time to datetime format
        fileDateTime = datetime([dateStr, timeStr], 'InputFormat', 'yyyyMMddHHmmss');

        % Check if the date already exists in the structure
        subjectID = findSubjectID(subjects, dateStr);

        % Create a new subject ID if necessary
        if isempty(subjectID)
            subjectID = ['s', num2str(subjectCounter)];
            subjectCounter = subjectCounter + 1;
        end

        % Append the current file data to the subject's file list
        if ~isfield(subjects, subjectID)
            subjects.(subjectID).files = struct('fileName', fileName, 'dateTime', fileDateTime);
        else
            subjects.(subjectID).files(end+1) = struct('fileName', fileName, 'dateTime', fileDateTime);
        end
    end

    % Determine the "ON" and "OFF" status for each subject based on time
    subjectNames = fieldnames(subjects);
    for i = 1:length(subjectNames)
        % Get the list of files and their datetimes
        subjectFiles = subjects.(subjectNames{i}).files;
        allTimes = [subjectFiles.dateTime];

        % Find the earliest and latest times
        [~, earliestIdx] = min(allTimes);
        [~, latestIdx] = max(allTimes);

        % Assign "ON" to the latest time and "OFF" to the earliest time
        subjects.(subjectNames{i}).ON = subjectFiles(latestIdx).fileName;
        subjects.(subjectNames{i}).OFF = subjectFiles(earliestIdx).fileName;

        % Remove the intermediate 'files' field
        subjects.(subjectNames{i}) = rmfield(subjects.(subjectNames{i}), 'files');
    end
end

% Helper function to find the subject ID based on the date
function subjectID = findSubjectID(subjects, dateStr)
    subjectID = '';  % Initialize subjectID as empty
    subjectNames = fieldnames(subjects);
    for i = 1:length(subjectNames)
        % Compare the date strings
        if strcmp(subjects.(subjectNames{i}).files(1).fileName(1:8), dateStr)
            subjectID = subjectNames{i};
            return;
        end
    end
end