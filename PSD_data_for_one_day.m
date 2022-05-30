%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                PSD based on data for one day                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save data as table and load table into MATLAB
clear; clc; close all;

%Save data to tables
% L = readtable('sid1005d2021-11-06.pkl.csv');
% T = readtable('sid1004d2021-11-06.pkl.csv');
% 
% save('transverse_and_longitudinal_data_2021-11-06','L','T')


%Load data for one day 
load('transverse_and_longitudinal_data_2021-11-06')

%% 

%Define the transverse and longitudinal sensor data
sensorT = T.value;   
sensorL = L.value;  

%Compute the Welch's PSD estimate
fs = 25;                                                                   % Sampling frequency (Hz)
N = 2^13;                                                                  % Number of points in frequency (blocksize)                                                               % Number of points in frequency (blocksize). FFT is desinged to be to a power of 2 (^2)
Nf = N/4;
win = hann(Nf);                                                            % Window type
nover = 3*Nf/4;                                                            % Number of overlapping samples

%Compute the Welch's PSD estimate
[SxxT,f1] = pwelch(sensorT,win,nover,Nf,fs);                               % Estimate PSD
[SxxL,f2] = pwelch(sensorL,win,nover,Nf,fs);                               % Estimate PSD

%Plot tranverse and longitidinal in same plot
figure()
plot(f1,10*log(SxxT)/(2*pi))
hold on 
plot(f2,10*log(SxxL)/(2*pi))
legend('Transverse','Longitudinal')
grid on 
title('Welch Power Spectral Density Estimate','Interpreter','latex','FontSize',13)
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',11)

%Plot tranverse and longitidinal in different plots
figure()
subplot(211)
plot(f1,10*log(SxxT))
title('Welch Power Spectral Density Estimate - Transverse')
xlabel('Frequency [Hz]')
grid on 
subplot(212)
plot(f2,10*log(SxxL))
title('Welch Power Spectral Density Estimate - Longitudinal')
xlabel('Frequency [Hz]')
grid on 
