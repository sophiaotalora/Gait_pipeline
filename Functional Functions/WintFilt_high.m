%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for implementing the forward and reverse 2nd order Butterworth digital 
% filter of % D.A. Winter as outlined in Winter, D.A., Biomechanics and Motor Control of Human 
% Movement, Wiley, 1990 pp36-43.
%
% Written by Dr. Gerard Lyons, June 20th 1998
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Filter_Output=WintFilt_high(Filter_Input, fcut_off, fsample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function requires as input an array contained the data to be filtered (Filter_Input), 
% the required cut-off frequency (fcut_off) and the sample frequency (fsample).
% The function outputs an array containing the filtered signal (Filter_Output).
% Format is: Filter_Output=WinFilt(Filter_Input, fcut_off, fsample)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program uses a cut-off freq of fcut_off in Hz, this value needs to be divided by
%  0.802 to get the effective cut-off due the effects of the forward and reverse passes 
%  through the filter		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc = fcut_off;%/0.802;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [B,A] = butter(2,Wn) designs an second order lowpass digital Butterworth filter and % returns the filter coefficients in length 3 vectors B and A. The cut-off frequency Wn
% must be 0.0 < Wn < 1.0, with 1.0 corresponding to half the sample rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wn=(fc*2)/(fsample);
[B,A] = butter(4,wn,'high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% FILTFILT Zero-phase forward and reverse digital filtering. Y = FILTFILT(B, A, X) 
% filters the data in vector X with the filter described by vectors A and B to create the  
% filtered data Y.  The filter is described by the difference equation:
% 
% 	  y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
% 	                   - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% 
% After filtering in the forward direction, the filtered sequence is then reversed and run 
% back through the filter. The resulting sequence has precisely zero-phase distortion and
% double the filter order. Care is taken to minimize startup and ending transients by  
% matching initial conditions. See also FILTER.			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Filter_Output = filtfilt(B, A, Filter_Input);
