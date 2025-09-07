clear all;
close all;
clc;
pkg load signal;
V_dc = 100;
R = 10;
L = 0.025;
f = 50;
T = 1/f;

dt = 10^(-6);

for m_f = [40 200] % Iterate for two values of m_f
    m_a = 0.9; % Amplitude modulation index
t = 0:dt:6*T-dt;

vacheck = m_a*sin(2*pi*f*t);
vbcheck = m_a*sin(2*pi*f*t - 2*pi/3);
vccheck = m_a*sin(2*pi*f*t + 2*pi/3);
Vtriangular = sawtooth(2*pi*f*m_f*t -pi/2,0.5);
q_1 = zeros(length(t), 1);
q_2 = zeros(length(t), 1);
q_3 = zeros(length(t), 1);
q_4 = zeros(length(t), 1);
q_5 = zeros(length(t), 1);
q_6 = zeros(length(t), 1);
Va = zeros(length(t), 1);
Vb = zeros(length(t), 1);
Vc = zeros(length(t), 1);
Vab = zeros(length(t), 1);
Vbc = zeros(length(t), 1);
Vca = zeros(length(t), 1);
Van = zeros(length(t), 1);
Vbn = zeros(length(t), 1);
Vcn = zeros(length(t), 1);
Ia = zeros(length(t),1);
Ib = zeros(length(t),1);
Ic = zeros(length(t),1);
for i = 1:1:length(t)
    if vacheck(i) > Vtriangular(i)
        q_1(i) = 1;
        Va(i) = V_dc;
    elseif vacheck(i) < Vtriangular(i)
        q_4(i) = 1;
        Va(i) = 0;
    end
    if vbcheck(i) > Vtriangular(i)
        q_3(i) = 1;
        Vb(i) = V_dc;
    elseif vbcheck(i) < Vtriangular(i)
        q_6(i) = 1;
        Vb(i) = 0;
    end
    if vccheck(i) > Vtriangular(i)
        q_2(i) = 1;
        Vc(i) = V_dc;
    elseif vccheck(i) < Vtriangular(i)
        q_5(i) = 1;
        Vc(i) = 0;
    end
    Vab(i) = Va(i) - Vb(i);
    Vbc(i) = Vb(i) - Vc(i);
    Vca(i) = Vc(i) - Va(i);
    Van(i) = (Vab(i)-Vca(i))/3;
    Vbn(i) = (Vbc(i)-Vab(i))/3;
    Vcn(i) = (Vca(i)-Vbc(i))/3;
    if i > 1
        Ia(i)=L/(R*dt+L)*Ia(i-1) + dt/(R*dt+L)*Van(i);
        Ib(i)=L/(R*dt+L)*Ib(i-1) + dt/(R*dt+L)*Vbn(i);
        Ic(i)=L/(R*dt+L)*Ic(i-1) + dt/(R*dt+L)*Vcn(i);
        I1(i) = Ia(i);
        I2(i) = Ib(i);
        I3(i) = Ic(i);
        q_1(i) = q_1(i)*I1(i);
        q_4(i) = q_4(i)*I1(i);
        q_3(i) = q_3(i)*I2(i);
        q_6(i) = q_6(i)*I2(i);
        q_2(i) = q_2(i)*I3(i);
        q_5(i) = q_5(i)*I3(i);

    end

end

% Signal names
figure
signals = {'Va', 'Vb', 'Vc', 'Vab', 'Vbc', 'Vca', 'Van', 'Vbn', 'Vcn'};
for i = 1:length(signals)
    subplot(3, 3, i)
    plot(t, eval(signals{i}))
    ylim([-110 110])
    xlim([-0.001 0.08])
    ylabel('V(V)')
    xlabel('t(sec)')
    title([signals{i}, ' for case mf =', num2str(m_f)])
end
% Signal names
figure;
signals = {'Ia', 'Ib', 'Ic'};
for i = 1:length(signals)
    subplot(3, 1, i)
    plot(t, eval(signals{i}))
    ylabel('I(A)')
    xlabel('t(sec)')
    title([signals{i}, ' for case mf = ', num2str(m_f)])
end

figure
% Subplot 1
subplot(3, 1, 1)
plot(t, q_1, 'LineWidth', 1.5)
hold on
plot(t, q_4, 'LineWidth', 1.5)
hold off
ylabel('I(A)')
xlabel('t(sec)')
legend('I q_1', 'I q_4')
title('Current in q_1 and q_4')
% Subplot 2
subplot(3, 1, 2)
plot(t, q_3, 'LineWidth', 1.5)
hold on
plot(t, q_6, 'LineWidth', 1.5)
hold off
ylabel('I(A)')
xlabel('t(sec)')
legend('I q_3', 'I q_6')
title('Current in q_3 and q_6')
% Subplot 3
subplot(3, 1, 3)
plot(t, q_2, 'LineWidth', 1.5)
hold on
plot(t, q_5, 'LineWidth', 1.5)
hold off
ylabel('I(A)')
xlabel('t(sec)')
legend('I q_2', 'I q_5')
title('Current in q_2 and q_5')
% Calculate RMS values
Van_rms = rms(Van(1:20000));
Ia_rms = rms(Ia(1:20000));
Vbn_rms = rms(Vbn(1:20000));
Ib_rms = rms(Ib(1:20000));
Vcn_rms = rms(Vcn(1:20000));
Ic_rms = rms(Ic(1:20000));
% Calculate Z, P_out, and P_in
poutput = Ia .* Van + Ib .* Vbn + Ic .* Vcn;
S = Van_rms * Ia_rms; % Any phase can be used
P_in = Ia_rms.^2 * R; % Any phase can be used
idc = sqrt(poutput / (sqrt(R^2 + (2*pi*f*L)^2)));
figure;
subplot(2,1,1)
plot(t, idc,'g', 'LineWidth', 1.5);
ylabel('I (A)');
title('Current dc');
xlabel('t (sec)');
subplot(2,1,2)
plot(t, poutput,'black', 'LineWidth', 1.5);
title('Power output ');
ylabel('P (W)');
xlabel('t (sec)');
% Add a horizontal line at y = 390.21
yValue = 66;
xLimits = xlim;  % Get the x-axis limits
line(xLimits, [yValue, yValue], 'Color', 'blue', 'LineStyle', '--');
%% power factor
powerfactor = P_in(1)/S
%%All the harmonics
NFFT = 20000 ;
frequency = (1/dt)/2*linspace(0,1,NFFT/2+1);
harmonics_signals = {Van, Vbn, Vcn, Vab, Vbc, Vca};
harmonics_titles = {'Van', 'Vbn', 'Vcn', 'Vab', 'Vbc', 'Vca'};
figure;
for i = 1:numel(harmonics_signals)
    signal = harmonics_signals{i};
    title_text = ['Harmonics for ' harmonics_titles{i}];
    Vharm = fft(signal, NFFT) / 20000;
    subplot(2, 3, i);
    plot(frequency, 2*abs(Vharm(1:NFFT/2+1)));
    title(title_text);
    xlabel('F(Hz)');
end

end
