

function y=filtering(x)
F_nyq=122/2;
[b,a]=butter(4,0.05/F_nyq,'high');%drift suppression (butterworth-high pass filter) to remove slow signal changes
x_filtrato=filtfilt(b,a,x);
[b,a]=butter(4,40/F_nyq,'low');%
x_filtrato2=filtfilt(b,a,x_filtrato);
y=x_filtrato2;



% load('M248.mat');
% x=yy;
% load('yy.mat');
%figure
%plot(x)
%title('Filtered signal')
%xlabel('time(msec)')
%ylabel('Voltage(mV)')
%% drift suppression (butterworth-high pass filter) to remove slow signal changes
%(cutoff frequency=1Hz)
%hold on
%plot(x_filtrato,'g')
%% butterworth filter (low pass) to eliminate frequencies higher than 30Hz
%(cutoff frequency=30Hz)
%hold on
%plot(x_filtrato2,'m')


