function [DominantFreq,AmplitudeNorm,WidthNorm,SlopeNorm,RangeNorm,IntegratedPower] = ...
    extractFreqBasedFeautres_AnerWeiss_norm(Data,Nfft,Fs,HanningWindow_Size,Settings,PlotSettings)

MaxRangePSD = Settings.MaxRangeFrequencyPSD;
MinRangePSD = Settings.MinRangeFrequencyPSD;

%% Measures tested by Weiss et al. 2013
% Spectral measures
for i=1:size(Data,2)
    % Band-filtering of accelerometry prior to gait features extraction
    %     Data(:,i) = LPFilter(Data(:,i),Fs,3,2); %15
    %     Data(:,i) = HPFilter(Data(:,i),Fs,0.5,2); %0.25
    AccWin_i = (Data(:,i)-mean(Data(:,i)))/std(Data(:,i));                  % normalize window
    w = hanning(floor(HanningWindow_Size));
    [PW_i,FW] = pwelch(AccWin_i(1:HanningWindow_Size),w,50,Nfft,Fs);        % normalized (window)    
    % Normalize PSD
    IntPW_i = cumsum(PW_i);       
    PW_i = PW_i/IntPW_i(end); %PW_i = PW_i/max(PW_i);
    if PlotSettings.SpectralFatures
        figure,plot(FW,PW_i)
    end
    IXFRange = find(FW>=MinRangePSD & FW<= MaxRangePSD);                                     % increased frequency range since ML needs a broader one
    FDind = IXFRange(find(PW_i(IXFRange)==max(PW_i(IXFRange)),1,'first'));
    FD = FW(FDind);
    FDAmp = PW_i(FDind);
    FDindRange = [find(PW_i<0.5*FDAmp & FW<FW(FDind),1,'last'), find(PW_i<0.5*FDAmp & FW>FW(FDind),1,'first')];
    if numel(FDindRange) == 2
        FDWidth = diff(FW(FDindRange));
    else
        FDWidth = nan;
    end
    if FDind ~= min(IXFRange) && FDind ~= max(IXFRange)
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
    AmplitudeNorm(1,i) = FDAmp;
    WidthNorm(1,i) = FDWidth;
    IntegratedPower(1,i) = IntPW_i(end); % Mannini
end
SlopeNorm = AmplitudeNorm./WidthNorm;
% Temporal measure
RangeNorm = max(Data)-min(Data);