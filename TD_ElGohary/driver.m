function [status] = driver(indir, outdir)
warning off
rmpath(genpath('./Library/'))
warning on
addpath(genpath('./Library/'))

% load the data
load(fullfile(indir,'data.mat'));

% load the GSD output
load(fullfile(indir,'GSDA_Output.mat'));

% select the sensor unit to use. All options can be found in
% ./Library/process_data.m
sensor_string = getenv_string('SENSOR', 'SU');
disp(strcat('Sensor: ', sensor_string));


% determine if the data set is from Free-living or Laboratory
% get the TimeMeasures in the data
time_measure_list = fieldnames(Data);
time_measure_list=GetFields(Data,'data',{});

% % look at the first time measure
% time_measure_i = Data.(time_measure_list{1});
%     
% % get a list of recordings
% recording_list = fieldnames(time_measure_i);

% check if there are "Recordings" or "Tests", get the fieldnames to iterate
% over accordingly
% if all(contains(recording_list, 'Recording'))
%     disp('Free living data set')
%     field_names = fieldnamesr(data,2);
% elseif all(contains(recording_list, 'Trial'))
%     disp('Lab data set')
%     field_names = fieldnamesr(data,3);
% else
%     error('Data set not recognized as free-living or lab data set')
% end
% Iterar sobre los sujetos
subject_names = fieldnames(Data);
for s = 1:length(subject_names)
    subject = subject_names{s};
    
    % Iterar sobre las condiciones (F, N)
    condition_names = fieldnames(Data.(subject));
    for c = 1:length(condition_names)
        condition = condition_names{c};
        
        % Iterar sobre los trials (Trial01, Trial02, ...)
        trial_names = fieldnames(Data.(subject).(condition));
        for t = 1:length(trial_names)
            trial = trial_names{t};
            
            % Extraer los datos del trial
            trial_data = Data.(subject).(condition).(trial);
            
            % Crear el campo SU combinando Acc, Gyr, Mag, y Fs
            if isfield(trial_data, 'Acc') && isfield(trial_data, 'Gyr') && ...
                    isfield(trial_data, 'Mag') && isfield(trial_data, 'Fs')
                
                Data.(subject).(condition).(trial).SU.LowerBack = struct( ...
                    'Acc', trial_data.Acc, ...
                    'Gyr', trial_data.Gyr, ...
                    'Mag', trial_data.Mag, ...
                    'Fs', trial_data.Fs);
            else
                % Crear un campo SU vacío si faltan datos
                Data.(subject).(condition).(trial).SU.LowerBack = struct('Acc', [], 'Gyr', [], 'Mag', [], 'Fs', []);
            end
            
            % Eliminar los campos originales (Acc, Gyr, Mag, Fs)
            Data.(subject).(condition).(trial) = rmfield(Data.(subject).(condition).(trial), {'Acc', 'Gyr', 'Mag', 'Fs'});
        end
    end
end
field_names = fieldnamesr(Data,4);
% call one common function for free-living and lab data
[output_struct, output_struct_json] = process_data(Data, GSD_Output, field_names, sensor_string);
%% save results
base_filename = strcat('TD_ElGohary', '_', sensor_string);

% .json file
output_struct_json = struct('TD_Output',output_struct_json);
json_string = jsonencode(output_struct_json);

% save json with regular name
filename = strcat(base_filename, '.json');
fid = fopen(fullfile(outdir,filename),'wt');
fprintf(fid,json_string);
fclose(fid);

% save json with generic name
fid = fopen(fullfile(outdir,'TD_Output.json'),'wt');
fprintf(fid,json_string);
fclose(fid);

% .mat file
TD_Output = output_struct;
filename = strcat(base_filename, '.mat');

% regular name
save(fullfile(outdir,filename),'TD_Output');

% generic name
save(fullfile(outdir,'TD_Output.mat'),'TD_Output');
status = 'ok';

% Metadata
add_metadata(outdir, 'metadata.txt', base_filename);
end