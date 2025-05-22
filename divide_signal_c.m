function Data =divide_signal_c(Data, TD_Output)
States = fieldnamesr(Data, 2);
States = unique(regexprep(States, '.*\.', ''));
SubjectNames = fieldnames(Data);
Trials= fieldnamesr(Data, 3);
Trials = unique(regexprep(Trials, '.*\.', ''));
for iSubject = 1:length(SubjectNames)
    for iState = 1:length(States)
        for iTrial = 1:length(Trials)
       
            Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB = struct('turnings', [], 'straightWalk', []);
            Turns = [TD_Output.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).SU.LowerBack.TD.Turn_Start];
            TurnsDur = [TD_Output.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).SU.LowerBack.TD.Turn_Duration];
            Fs = Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).Fs.Acc; 
            total_samples = length(Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).Acc); 
            start_index = 1;                                               % Index start
    
            TurnsSize = length(Turns);
            for iTurn = 1:TurnsSize
               
                end_index = round(Turns(iTurn) * Fs);                      % Index end of straight gait before turn
    
                                                                           % Add straight gait interval to strc
                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.straightWalk = [ ...
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.straightWalk; ...
                    struct('start', start_index, 'end', end_index) ...
                ];
    
                %Index start and end of turning
                turning_start = end_index + 1;
                turning_end_index = min(turning_start + round(TurnsDur(iTurn) * Fs), total_samples);
    
                %Add turning if it's valid
                if turning_end_index > turning_start
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.turnings = [ ...
                        Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.turnings; ...
                        struct('start', turning_start, 'end', turning_end_index) ...
                    ];
                end
    
                start_index = turning_end_index + 1;                       % Update the start for next straight walk
            end
    
            
            if start_index <= total_samples                                % Last straight walk before turn
                Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.straightWalk = [ ...
                    Data.(SubjectNames{iSubject}).(States{iState}).(Trials{iTrial}).WB.straightWalk; ...
                    struct('start', start_index, 'end', total_samples) ...
                ];
            end
        end
    end
end
