function subjects = processfilesh5_2(DirectoryPath_Data)
    folderPath = DirectoryPath_Data;
    
    files = dir(fullfile(folderPath, '*.h5'));
    
    subjects = struct();
    
    % Iterate through the file pairs
    numSubjects = length(files) / 2;
    
    for i = 1:numSubjects
        offIndex = (i - 1) * 2 + 1;                                        % Calculate indices for OFF and ON files
        onIndex = offIndex + 1;
        
        subjectName = ['s' num2str(i)];                                    % Create subject identifier, e.g., 's1', 's2', ...  
     
        subjects.(subjectName).Test_OFF = files(offIndex).name;            % Assign files to OFF and ON fields
        subjects.(subjectName).Test_ON = files(onIndex).name;
    end

end