clear all; % Clear all variables in the workspace
close all; % Close all open figures
clc; % Clear the command window

pkg load signal; % Load the signal package for additional signal processing functions

V_dc = 100; % DC voltage
R = 10; % Resistance
L = 0.025; % Inductance
f = 50; % Frequency
T = 1/f; % Period

dt = 10^(-6); % Time step

for m_f = [40 200] % Iterate for two values of m_f
    m_a = 0.9; % Amplitude modulation index
    t = 0:dt:6*T-dt; % Time vector

    Vsine = m_a*sin(2*pi*f*t); % Sine waveform with amplitude modulation
    Vsineopp = -Vsine; % Opposite of the sine waveform

    Vsaw = sawtooth(2*pi*f*m_f*t - pi/2, 0.5); % Sawtooth waveform with frequency modulation

    q_1 = zeros(length(t), 1); % Initialize q1
    q_2 = zeros(length(t), 1); % Initialize q2
    q_3 = zeros(length(t), 1); % Initialize q3
    q_4 = zeros(length(t), 1); % Initialize q4

    Va = zeros(length(t), 1); % Initialize Va
    Vb = zeros(length(t), 1); % Initialize Vb
    Voutput = zeros(length(t), 1); % Initialize the output voltage
    Ioutput = zeros(length(t), 1); % Initialize the output current

    for i = 1:1:length(t)
        % Determine the states of q1, q2, q3, and q4 based on the comparison between Vsine and Vsaw
        if Vsine(i) > Vsaw(i)
            q_1(i) = 1;
            Va(i) = V_dc;
        elseif Vsine(i) < Vsaw(i)
            q_4(i) = 1;
            Va(i) = 0;
        end

        % Determine the states of q2 and q3 based on the comparison between Vsineopp and Vsaw
        if Vsineopp(i) > Vsaw(i)
            q_3(i) = 1;
            Vb(i) = V_dc;
        elseif Vsineopp(i) < Vsaw(i)
            q_2(i) = 1;
            Vb(i) = 0;
        end

        Voutput(i) = Va(i) - Vb(i); % Calculate the output voltage

        if i > 1
            % Calculate the output current using the previous current and the output voltage
            Ioutput(i) = L/(R*dt+L)*Ioutput(i-1) + dt/(R*dt+L)*Voutput(i);
            I(i) = Ioutput(i);

            % Multiply the states q1, q2, q3, and q4 by the current
            q_1(i) = q_1(i) * I(i);
            q_2(i) = q_2(i) * I(i);
            q_3(i) = q_3(i) * I(i);
            q_4(i) = q_4(i) * I(i);
        end
    end

figure;
% Plot Voutput
subplot(2,1,1);
plot(t, Voutput, 'c', 'LineWidth', 1.5);
title(['Vout for mf = ' num2str(m_f)]);
ylabel('V(V)');
xlabel('t(sec)');
xlim([-0 0.08]);
ylim([-130 130]);


subplot(2,1,2)
plot(t,Ioutput, 'g', 'LineWidth', 2);
ylabel('I(A)');
xlabel('t(sec)');
title(['Iout for mf = ' num2str(m_f)]);


figure
subplot(2,1,1)
plot(t,q_1, 'm', 'LineWidth', 1.5);
hold on;
plot(t,q_4, 'g--', 'LineWidth', 1.5);
hold off;
legend('I q_1','I q_4');
title(['Iq_1 and Iq_4 for mf = ' num2str(m_f)]);
ylabel('I(A)');
xlabel('t(sec)');
subplot(2,1,2)
plot(t,q_2, 'b', 'LineWidth', 1.5);
hold on;
plot(t,q_3, 'r--', 'LineWidth', 1.5);
hold off;
legend('I q_2','I q_3');
title(['Iq_2 and Iq_3 for mf = ' num2str(m_f)]);
ylabel('I(A)');
xlabel('t(sec)');

Vout_rms = rms(Voutput(1:20000));
Iout_rms = rms(Ioutput(1:20000));

Z = sqrt(R^2 + (2*pi*f*L)^2);
poutput = Ioutput .* Voutput;
S = Vout_rms * Iout_rms;
P = Iout_rms^2 * R * ones(20000, 1);
idc = sqrt(poutput / Z);

figure;
subplot(2,1,1)
plot(t, idc,'g', 'LineWidth', 1.5);
ylabel('I (A)');
title(['I_{DC} for m_f = ' num2str(m_f)]);
xlabel('t (sec)');
subplot(2,1,2)
plot(t, poutput,'black', 'LineWidth', 1.5);
title(['Power output for m_f = ' num2str(m_f)]);
ylabel('P (W)');
xlabel('t (sec)');
% Add a horizontal line at y = 390.21
yValue = 262.63;
xLimits = xlim;  % Get the x-axis limits
line(xLimits, [yValue, yValue], 'Color', 'blue', 'LineStyle', '--');

%% power factor
pf = P / S;
disp(['PowerFactor for m_f = ' num2str(m_f) ': ']);
disp(pf(1));


NFFT = 20000;
frequency = (1/ dt) / 2 * linspace(0, 1, NFFT/2+1);
Vharm = fft(Voutput, NFFT) / 20000;

figure;
plot(frequency, 2 * abs(Vharm(1:NFFT/2+1)));
title(['Harmonics for m_f = ' num2str(m_f)]);
xlabel('F (Hz)');
xlim([0 200000]);


end
