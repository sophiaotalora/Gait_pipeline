function [DominantFreq,Amplitude,Width,Slope,Range,MeanPower,MedianPower,MeanFreq,MedianFreq] = ...
    extractFreqBasedFeautres_AnerWeiss(Data,Nfft,Fs,HanningWindow_Size,Settings,PlotSettings)

MaxRangePSD = Settings.MaxRangeFrequencyPSD;
MinRangePSD = Settings.MinRangeFrequencyPSD;

%% Measures tested by Weiss et al. 2013
% Spectral measures
for i=1:size(Data,2)
    % Band-filtering of accelerometry prior to gait features extraction
    %     Data(:,i) = LPFilter(Data(:,i),Fs,3,2); %15
    %     Data(:,i) = HPFilter(Data(:,i),Fs,0.5,2); %0.25
    % Calculation
    AccWin_i = (Data(:,i)-mean(Data(:,i)))/std(Data(:,i));                  % normalize window
    w = hanning(floor(HanningWindow_Size));
    [PW_i,FW] = pwelch(AccWin_i(1:HanningWindow_Size),w,50,Nfft,Fs);        %normalized (window)
    IXFRange = find(FW>=MinRangePSD & FW<= MaxRangePSD); % frequency range, between 0.3 and 10 Hz % increased frequency range since ML needs a broader one
    FDind = IXFRange(find(PW_i(IXFRange)==max(PW_i(IXFRange)),1,'first'));  % sample at which the peak of Power is found
    FD = FW(FDind); % frequency of the peak PSD
    FDAmp = PW_i(FDind); % amplitude of the peak PSD
    FDindRange = [find(PW_i<0.5*FDAmp & FW<FW(FDind),1,'last'), find(PW_i<0.5*FDAmp & FW>FW(FDind),1,'first')];
    if PlotSettings.SpectralFatures
        figure,plot(FW,PW_i)
    end
    % range of samples, half amplitude before and after peak
    if numel(FDindRange) == 2
        FDWidth = diff(FW(FDindRange));
    else
        FDWidth = nan;
    end
    if FDind ~= min(IXFRange) && FDind ~= max(IXFRange) % if  sample of peak is not the minimum or maximum of the predefined frequency range
        %% ParabolaVertex
        VertexIX = [-1 0 1] + FDind;
        [FD,FDAmp] = ParabolaVertex(FW(VertexIX),PW_i(VertexIX));
        %%
        FDindRange = [find(PW_i<0.5*FDAmp & FW<FD,1,'last'), find(PW_i<0.5*FDAmp & FW>FD,1,'first')];
        if numel(FDindRange) == 2
            StartP = PW_i(FDindRange(1)+[0 1]);
            StartF = FW(FDindRange(1)+[0 1]);
            StopP = PW_i(FDindRange(2)-[0 1]);
            StopF = FW(FDindRange(2)-[0 1]);
            FDRange = [interp1(StartP,StartF,0.5*FDAmp) , interp1(StopP,StopF,0.5*FDAmp)];
            FDWidth = diff(FDRange);
        end
    end
    DominantFreq(1,i) = FD;
    Amplitude(1,i) = FDAmp;
    Width(1,i) = FDWidth;
    MeanPower(1,i) = mean(PW_i); % Qu Wei et al.
    MedianPower(1,i) = median(PW_i); % Qu Wei et al.
    MeanFreq(1,i) = mean(FW); % Qu Wei et al.
    MedianFreq(1,i) = median(FW); % Qu Wei et al.
end
Slope = Amplitude./Width;
% Temporal measure
Range = max(Data)-min(Data);

