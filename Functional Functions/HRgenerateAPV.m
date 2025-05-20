function [  HRAPV ] = HRgenerateAPV( f1a )
%Generate the harmonic ratio from the accelleration data for the AP and the
%V directions
%

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

% Calulation of the Harmonic ratio
% AP and V compnents = 
% sum of Amplitudes of even harmonics/sum of Amplitudes of odd harmonics

HRAPV = sum(FourierResult(2:2:end))/sum(FourierResult(1:2:end));




