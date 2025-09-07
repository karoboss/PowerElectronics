%Starting commands
close all;
clear all;
clc;
pkg load signal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   2.FFT PLOTS AND BASIC PLOTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = 10^(-5); % Sample time
Fs = 1/Ts; % Sampling frequency
L = 20000; % Length of signal

%%Load functions for the Signals to apply Fourier transformation
t = load('time.mat');

y = load('voltage1.mat');

y2 = load('voltage2.mat');

y3 = load('voltage3.mat');

y_curr = load('current.mat');

l4=length(y.voltage1);
%Plot our basic signals

% Plot Voltage 1
figure
plot(t.time,y.voltage1)
title('Voltage 1')
xlabel('time (sec)')
ylabel('Voltage 1')

% Plot Voltage 2
figure
plot(t.time,y2.voltage2)
title('Voltage 2')
xlabel('time (sec)')
ylabel('Voltage 2')

% Plot Voltage 3
figure
plot(t.time,y3.voltage3)
title('Voltage 3')
xlabel('time (sec)')
ylabel('Voltage 3')

% Plot Current
figure
plot(t.time,y_curr.current)
title('Current')
xlabel('time (sec)')
ylabel('Current')

%Fourier our basic signals
NFFT = 20000; %nfft sampling
f = Fs/2 * linspace(0,1,NFFT/2+1);%Frequency vector


% Fourier transformation Voltage 1
Y = fft(y.voltage1,NFFT)/L;
% Fourier transformation Voltage 2
Y2 = fft(y2.voltage2,NFFT)/L;
% Fourier transformation Voltage 3
Y3 = fft(y3.voltage3,NFFT)/L;
% Fourier transformation Current
Y_CURR = fft(y_curr.current,NFFT)/L;

% Plot Fourier Voltage 1
figure
plot(f(1:250), abs(Y(1:250)));
title('Fourier Transform Voltage 1')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% Plot Fourier Voltage 2
figure
plot(f(1:250), abs(Y2(1:250)));
title('Fourier Transform Voltage 2')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% Plot Fourier Voltage 3
figure
plot(f(1:250), abs(Y3(1:250)));
title('Fourier Transform Voltage 3')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

% Plot Fourier Current
figure
plot(f(1:250), abs(Y_CURR(1:250)));
title('Fourier Transform Current')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   3.SPECTRUM PLOTS     %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot single-sided amplitude spectrum voltage 1.
figure
plot(f(1:250),2 * abs(Y(1:NFFT/2+1))(1:250))
title('Single-Sided Amplitude Spectrum of Voltage 1 with NFFT = 20000')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%Plot single-sided amplitude spectrum voltage 2.
figure
plot(f(1:250),2 * abs(Y2(1:NFFT/2+1))(1:250))
title('Single-Sided Amplitude Spectrum of Voltage 2 with NFFT = 20000')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%Plot single-sided amplitude spectrum voltage 3.
figure
plot(f(1:250),2 * abs(Y3(1:NFFT/2+1))(1:250))
title('Single-Sided Amplitude Spectrum of Voltage 3 with NFFT = 20000')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%Plot single-sided amplitude spectrum current.
figure
plot(f(1:250),2 * abs(Y_CURR(1:NFFT/2+1))(1:250))
title('Single-Sided Amplitude Spectrum of Current with NFFT = 20000')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%%%%%%%%%%%%%%%%%%%%%%%%%%% 4. THD AND OTHER CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  VOLTAGE1 THD  %%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:1/Fs:1-1/Fs; % Time vector (seconds)
mag_voltage1 = abs(Y);

%Here we can calculate the fundamental frequency at the first harmonic
[max_mag1, fund_index1] = max(mag_voltage1);

%Here we are using a for loop iterating over the harmonics and putting in a variable both the sums of the voltages in power of 2
%and in power of 1
total_harmonic_distortion1 = 0;
dpf1=0;
for n = 2:length(mag_voltage1)/2 % Iterate over harmonic components
    total_harmonic_distortion1 = total_harmonic_distortion1 + mag_voltage1(n)^2; % Add square of magnitude of harmonic component
    dpf1=dpf1+mag_voltage1(n);
end
dpf1=dpf1;
%Here we use the type as shown in slides 3 of presentation to calculate the thd
total_harmonic_distortion1 =100* sqrt(total_harmonic_distortion1-max_mag1.^2)/max_mag1;% thd type
%We use the same method for the rest of the signals

%%%%%%%%%%%%%%%%%%%%%%%   VOLTAGE2 THD   %%%%%%%%%%%%%%%%%%%%%%%%%%%
mag_voltage2 = abs(Y2);
%%line to find the fundamental freq
[max_mag2, fund_index2] = max(mag_voltage2);
% calculate thd using for lopp
total_harmonic_distortion2 = 0;
dpf2=0;
for n = 2:length(mag_voltage2)/2 % Iterate over harmonic components
    total_harmonic_distortion2 = total_harmonic_distortion2 + mag_voltage2(n)^2; % Add square of magnitude of harmonic component
    dpf2=dpf2+mag_voltage2(n);
end
total_harmonic_distortion2 = 100*sqrt(total_harmonic_distortion2-max_mag2.^2)/max_mag2; % thd type

%%%%%%%%%%%%%%%%%%%%%%%  VOLTAGE3 THD  %%%%%%%%%%%%%%%%%%%%%%%%%%%

mag_voltage3 = abs(Y3);
%%line to find the fundamental freq
[max_mag3, fund_index3] = max(mag_voltage3);
% calculate thd using for lopp
total_harmonic_distortion3 = 0;
dpf3=0;
for n = 2:length(mag_voltage3)/2 % Iterate over harmonic components
    total_harmonic_distortion3 = total_harmonic_distortion3 + mag_voltage3(n)^2; % Add square of magnitude of harmonic component
    dpf3=dpf3+mag_voltage3(n);
end
total_harmonic_distortion3 =100* sqrt(total_harmonic_distortion3-max_mag3.^2)/max_mag3; % thd type

%%%%%%%%%%%%%%%%%%%%%%%  CURRENT THD  %%%%%%%%%%%%%%%%%%%%%%%%%%%
mag_voltage4 = abs(Y_CURR);
%%line to find the fundamental freq
[max_mag4, fund_index4] = max(mag_voltage4);
% calculate thd using for lopp
total_harmonic_distortion4 = 0;
dpf4=0;
for n = 2:length(mag_voltage4)/2 % Iterate over harmonic components
    total_harmonic_distortion4 = total_harmonic_distortion4 + mag_voltage4(n)^2; % Add square of magnitude of harmonic component
    dpf4=dpf4+mag_voltage4(n);
end
total_harmonic_distortion4 =100* sqrt(total_harmonic_distortion4-max_mag4.^2)/max_mag4; % thd type

%%%%%%%%%%%%%%%%%%%%%%%   ALL DISPLAYS   %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display thd
disp(['%%%%%%%%%%%%%%%%%%%%%%%  THD CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);

disp(['The value of total harmonic distortion (thd) for voltage1 is:' num2str(total_harmonic_distortion1)]);
disp(['The value of total harmonic distortion (thd) for voltage2 is:' num2str(total_harmonic_distortion2)]);
disp(['The value of total harmonic distortion (thd) for voltage3 is:' num2str(total_harmonic_distortion3)]);
disp(['The value of total harmonic distortion (thd) for current is:' num2str(total_harmonic_distortion4)]);


%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE DPF  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here we use a function to calculate the index at the fundamental harmonic frequency
%and then we use it to calculate the dpf which is the cos of the difference in angle
%Between every voltage signal and the current we have

[~, idx] = min(abs(f - 10));

disp(['%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE DPF  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);
disp(['The position of ', num2str(10), ' Hz in the frequency vector is ', num2str(idx)]);
DPF1 =  cos(angle(Y(idx))-(angle(Y_CURR(idx))));
DPF2 =  cos(angle(Y2(idx))-(angle(Y_CURR(idx))));
DPF3 =  cos(angle(Y3(idx))-(angle(Y_CURR(idx))));

%%%%%%%%%%%%%%%%%%%%%%%  ALL DISPLAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['The DPF for voltage1 is:' num2str(DPF1)]);
disp(['The DPF for voltage2 is:' num2str(DPF2)]);
disp(['The DPF for voltage3 is:' num2str(DPF3)]);

%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE  POWER FACTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we use the type to calculate the power factor based on the
%type Ibaseharmonic/Isharmonic*dpf of each voltage


%Calculation of the rms of the current
mag5=abs(Y_CURR);
counter5=0;
for n = 1:NFFT/2 +1
    counter5 = counter5 + 2*(mag5(n)^2);

end
rmscurr2=sqrt(sum(counter5));

%Actual type for power factor
PF1 = (max_mag4/ rmscurr2) * DPF1;
PF2 = (max_mag4/ rmscurr2) * DPF2;
PF3 = (max_mag4/ rmscurr2) * DPF3;

%%%%%%%%%%%%%%%%%%%%%%%  ALL DISPLAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE  POWER FACTOR  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);
disp(['The Power Factor for voltage1 is:' num2str(PF1)]);
disp(['The Power Factor for voltage2 is:' num2str(PF2)]);
disp(['The Power Factor for voltage3 is:' num2str(PF3)]);

%%%%%%%%%%%%%%%%%%%%%%%   CALCULATE FAINOMENES  %%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['%%%%%%%%%%%%%%%%%%%%%%%  APPARENT POWER CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);

%%%%%%%%%%%%%%%%%%%%%%%   VOLTAGE1 RMS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%rms calculation based on sum of magnitudes in the power of 2 and all ins square root

mag1=abs(Y);
counter1=0;
for n = 1:NFFT/2 +1
    counter1 = counter1 + 2*(mag1(n)^2);

end
rmsv1=sqrt(sum(counter1));
%%%%%%%%%%%%%%%%%%%%%%%   VOLTAGE2 RMS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

mag2=abs(Y2);
counter2=0;
for n = 1:NFFT/2 +1
    counter2 = counter2 + 2*(mag2(n)^2);

end
rmsv2=sqrt(sum(counter2));
%%%%%%%%%%%%%%%%%%%%%%%   VOLTAGE3 RMS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

mag3=abs(Y3);
counter3=0;
for n = 1:NFFT/2 +1
    counter3 = counter3 + 2*(mag3(n)^2);

end
rmsv3=sqrt(sum(counter3));
%%%%%%%%%%%%%%%%%%%%%%%   CURRENT RMS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

mag4=abs(Y_CURR);
counter4=0;
for n = 1:NFFT/2 +1
    counter4 = counter4 + 2*(mag4(n)^2);

end
rmscurr=sqrt(sum(counter4));

%%%%%%%%%%%%%%%%%%%%%%%  ALL DISPLAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['RMS value of voltage signal 1 = ' num2str(rmsv1) ' V']);
disp(['RMS value of voltage signal 2 = ' num2str(rmsv2) ' V']);
disp(['RMS value of voltage signal 3 = ' num2str(rmsv3) ' V']);
disp(['RMS value of current signal  = ' num2str(rmscurr) ' V']);
%%%%%%%%%%%%%%%%%%%%%%%   FINAL CALCULATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%
S1 = rmsv1 * rmscurr;
S2 = rmsv2 * rmscurr;
S3 = rmsv3 * rmscurr;

%%FOR EVALUATION
%S1 = rms(y.voltage1) * rms(y_curr.current);
%S2 = rms(y2.voltage2) * rms(y_curr.current);
%S3 = rms(y3.voltage3) * rms(y_curr.current);

%%%%%%%%%%%%%%%%%%%%%%%  ALL DISPLAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['The Apparent Power(S) for voltage1 is:' num2str(S1)]);
disp(['The Apparent Power(S) for voltage2 is:' num2str(S2)]);
disp(['The Apparent Power(S) for voltage3 is:' num2str(S3)]);

%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE REAL POWER  %%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE REAL POWER  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);
%real power calculation based on the type from slides 3
P1 = rmsv1 * max_mag4 * DPF1;
P2 = rmsv2 * max_mag4 * DPF2;
P3 = rmsv3 * max_mag4 * DPF3;

%%%%%%%%%%SOS%%%%%%%%%%%%%%%%%%%%%%
%%HERE I TRIED USING THE RMS CURRENT NOT THE FUNDAMENTAL RMS AND THE PF WERE
%%VERY CLOSE TO THE REAL VALUES BUT THE TYPE SAID WITH THE FUNDAMENTAL FREQUENCY SO
%% I LET THAT BE AND PUT THIS IN COMMENTS
%P1 = rmsv1 * rmscurr * DPF1;
%P2 = rmsv2 * rmscurr * DPF2;
%P3 = rmsv3 * rmscurr * DPF3;


%%%%%%%%%%%%%%%%%%%%%%%  ALL DISPLAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['The Real Power(P) for voltage1 is:' num2str(P1)]);
disp(['The Real Power(P) for voltage2 is:' num2str(P2)]);
disp(['The Real Power(P) for voltage3 is:' num2str(P3)]);

%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE ISNTANT POWER  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%I had to reload cause i had some errors im not showing because they are vectors of 20000 length
load('voltage1.mat');
load('voltage2.mat');
load('voltage3.mat');
load('current.mat');
P_inst_1 = voltage1 .* current;
P_inst_2 = voltage2 .* current;
P_inst_3 = voltage3 .* current;

%%%%%%%%%%%%%%%%%%%%%%%  ALL DISPLAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(['%%%%%%%%%%%%%%%%%%%%%%%  CALCULATE ISNTANT POWER  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);
%disp(['The Instant Power(P_inst) for voltage1 is:' num2str(P_inst_1)]);
%disp(['The Instant Power(P_inst) for voltage2 is:' num2str(P_inst_2)]);
%disp(['The Instant Power(P_inst) for voltage3 is:' num2str(P_inst_3)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  5    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RMS CALCULATIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['%%%%%%%%%%%%%%%%%%%%%%%  RMS CALCULATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);
%%RMS CALCULATION FROM DEFINITION WITH 2 DIFFERENT METHODS

%%FIRST ONE
V2_rms = sqrt(sum(y2.voltage2.^2)/length(y2.voltage2));
disp(['The rms value from definition is:(FIRST TYPE)' num2str(V2_rms) ' V']);

%%SECOND ONE ALREADY CALCULATED FROM BEFORE CHECK CODE
disp(['The rms value from definition is:(SECOND TYPE)' num2str(rmsv2) ' V']);

%%RMS FROM THE SPECTRUM OF THE DIAGRAMM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RMS UP UNTIL 700HZ(FIRST WAY)    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_cutoff = 700;  % Cutoff frequency in Hz
ind_cutoff = find(f <= f_cutoff, 1, 'last');
magtest=abs(Y2(1:ind_cutoff));
countertest=0;
for n = 1:ind_cutoff
    countertest = countertest + 2*(magtest(n)^2);

end
rms1stway =sqrt(sum(countertest));
disp(['Rms value for voltage2 from spectrum: until 700hZ (1st way) ' num2str(rms1stway) ' V']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RMS UP UNTIL 700HZ (SECOND WAY)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttest= 0:1/Fs:1-1/Fs;  % Time vector in seconds
Ntest = length(y2);
mag_test=abs(Y2);
ftest=(0:Ntest-1)*(Fs/Ntest);
V2_rms_diagram = sqrt(sum(2*mag_test(1:ind_cutoff).^2));
disp(['Rms value for voltage2 from spectrum: until 700hZ (2nd way) ' num2str(V2_rms_diagram) ' V']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



