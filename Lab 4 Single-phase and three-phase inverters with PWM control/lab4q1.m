% Clear all variables, figures, and command window
clear all;
close all;
clc;
% Load the signal processing package
pkg load signal;
% Define the given parameters
Vin = 100;
R = 10;
L = 0.025;
f = 50;
T = 1/f;
dt = 2*10^(-5);
% Loop over different values of 'a' (30, 18)
for a = 30:-12:18
t = 0:dt:6*T-dt;
% Initialize array
Vin = zeros(length(t), 1);
Voutput = zeros(length(t), 1);
Ioutput = zeros(length(t), 1);
% Q initialization
q_1 = zeros(length(t), 1);
q_2 = zeros(length(t), 1);
q_3 = zeros(length(t), 1);
q_4 = zeros(length(t), 1);
% D initialization
d_1 = zeros(length(t), 1);
d_2 = zeros(length(t), 1);
d_3 = zeros(length(t), 1);
d_4 = zeros(length(t), 1);

amod = mod(a, 360) / 360 * T;
half_period=T/2;

for i = 1:1:length(t)
checking_mod=mod(t(i),T);
Vin(i) = 100;
% Determine the state of switches and output voltage based on the time period
if checking_mod > 0 && checking_mod< amod
    d_4(i) = 1;
    q_2(i) = 1;
    Voutput(i) = 0;
elseif checking_mod>= amod && checking_mod< half_period - amod
    q_1(i) = 1;
    q_2(i) = 1;
    Voutput(i) = 1;
elseif checking_mod>= half_period - amod && checking_mod< half_period
    q_1(i) = 1;
    d_3(i) = 1;
    Voutput(i) = 0;
elseif checking_mod>= half_period &&checking_mod< half_period + amod
    d_1(i) = 1;
    q_3(i) = 1;
    Voutput(i) = 0;
elseif checking_mod>= half_period + amod && checking_mod< 2*half_period - amod
    q_3(i) = 1;
    q_4(i) = 1;
    Voutput(i) = -1;
elseif checking_mod>= 2*half_period - amod &&checking_mod< 2*half_period
    d_2(i) = 1;
    q_4(i) = 1;
    Voutput(i) = 0;
end
% Calculate output voltage
Voutput(i) = Voutput(i)*Vin(i);

    if i > 1
        % Calculate output current
        Ioutput(i)=L/(R*dt+L)*Ioutput(i-1) + dt/(R*dt+L)*Voutput(i);
        I(i) = Ioutput(i);

        d_1(i) = d_1(i)*I(i);
        d_2(i) = d_2(i)*I(i);
        d_3(i) = d_3(i)*I(i);
        d_4(i) = d_4(i)*I(i);
        q_1(i) = q_1(i)*I(i);
        q_2(i) = q_2(i)*I(i);
        q_3(i) = q_3(i)*I(i);
        q_4(i) = q_4(i)*I(i);

    end

end

figure;
% Plot Voutput
subplot(2,1,1);
plot(t, Voutput, 'c', 'LineWidth', 1.5);
hold on;
plot(t, Vin, 'y--', 'LineWidth', 1.5);
hold off;
legend('V ac', 'V in');
title(['Vout for a = ' num2str(a)]);
ylabel('V(V)');
xlabel('t(sec)');
xlim([-0 0.1]);
ylim([-130 130]);

subplot(2,1,2)
plot(t,Ioutput, 'g', 'LineWidth', 2);
ylabel('I(A)');
xlabel('t(sec)');
title(['Iout for a = ' num2str(a)]);

figure
subplot(2,1,1);
plot(t,q_1, 'm', 'LineWidth', 1.5);
hold on;
plot(t,q_4, 'g--', 'LineWidth', 1.5);
hold off;
legend('I q_1','I q_4');
title(['Iq_1 and Iq_4 for a = ' num2str(a)]);
ylabel('I(A)');
xlabel('t(sec)');
subplot(2,1,2);
plot(t,q_2, 'b', 'LineWidth', 1.5);
hold on;
plot(t,q_3, 'r--', 'LineWidth', 1.5);
hold off;
legend('I q_2','I q_3');
title(['Iq_2 and Iq_3 for a = ' num2str(a)]);
ylabel('I(A)');
xlabel('t(sec)');

figure
% Plot D1, D2, D3, and D4
plot(t, d_1, 'b', 'LineWidth', 1.5);
hold on;
plot(t, d_2, 'r--', 'LineWidth', 1.5);
plot(t, d_3, 'g-.', 'LineWidth', 1.5);
plot(t, d_4, 'm:', 'LineWidth', 1.5);
hold off;
legend('I d_1','I d_2','I d_3','I d_4');
title(['Id_1 Id_2 Id_3 and Id_4 for a = ' num2str(a)]);
ylabel('I(A)');
xlabel('t(sec)');

% Calculate RMS values, power, and power factor
Vout_rms = rms(Voutput(1:1000));
Iout_rms = rms(Ioutput(1:1000));
Z = sqrt(R^2 + (2*pi*f*L)^2);
poutput = Ioutput .* Voutput;
S = Vout_rms * Iout_rms;
P = Iout_rms^2 * R * ones(1000, 1);
idc = sqrt(poutput / Z);

figure;
subplot(2,1,1);
plot(t, idc,'g', 'LineWidth', 1.5);
ylabel('I (A)');
title(['I_{DC} for a = ' num2str(a)]);
xlabel('t (sec)');
subplot(2,1,2);
plot(t, poutput,'black', 'LineWidth', 1.5);
ylabel('P (W)');
title(['Power output for a = ' num2str(a)]);
xlabel('t (sec)');
% Add a horizontal line at y = 390.21
yValue = 390.21;
xLimits = xlim;  % Get the x-axis limits
line(xLimits, [yValue, yValue], 'Color', 'blue', 'LineStyle', '--');


% Power factor calculation
pf = P / S;
disp(['PowerFactor for a = ' num2str(a) ': ']);
disp(pf(1));

NFFT = 1000;
frequency = (1 / dt) / 2 * linspace(0, 1, NFFT/2+1);
V1 = fft(Voutput, NFFT) / 1000;

figure;
plot(frequency, 2 * abs(V1(1:NFFT/2+1)));
xlim([-100 1000]);
xlabel('F (Hz)');
title(['Harmonics for a = ' num2str(a)]);

end



