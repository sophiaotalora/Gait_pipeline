function [ HRML ] = HRgenerateML( f1a )
%Generate the harmonic ratio from the accelleration data for the ML direction

na = 20;

%Fourier analysis 
f1=fft(f1a)/length(f1a);

%mean
mean1=f1(1);

%harmonics (amplitude and phase)

for h=1:na,
    mod1(h)=(abs(f1(h+1))*2);     
end;

FourierResult = mod1;


% ML Compnent
% sum of Amplitudes of odd harmonics/ sum of Amplitudes of even harmonics

HRML = sum(FourierResult(1:2:end))/sum(FourierResult(2:2:end));


