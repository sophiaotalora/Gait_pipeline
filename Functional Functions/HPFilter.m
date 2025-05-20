% ---------------------------------------------------------------------------
%     File:       HPfilter
%     Author:     
% 
%     Date:       
% 
%     Revision:   
%       06092004.1 Applied GaitTest header, RS
%
%
%
%     Purpose:
%
% ---------------------------------------------------------------------------
function[smooth,res]=HPFilter(data,Fs,Fc, order)

if (nargin==3)
    order = 2;
end

% HPFILTER.M
% HPfilter(data,Fs,Fc)
% Function file used for high pass-filtering with zero phase-shift.
% A butterworth filter is used. Fs is sample frequency, Fc is cut-off frequency
% This function requires the signal processing toolbox which contains "filtfilt",
% "butter" and other tools.
%
% Wiebren Zijlstra, 26 Augustus 1999

[m,n]=size(data);
if m<=3*order
    smooth=data;
    res=sqrt((1/m)*sum((data-smooth).^2));
else
    Fs=Fs./2;
    %Fc=Fc/.802;

    [B,A] = butter(order,Fc/Fs,'high'); % 2nd order high pass filter.

    % 4th order high-pass filter, zero phase lag.
    for count=1:n
        smooth(:,count) = filtfilt(B,A,data(:,count));
    end
    
    res=sqrt((1/m)*sum((data-smooth).^2));
end


