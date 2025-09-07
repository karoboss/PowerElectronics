clear all;
close all;
clc;
pkg load signal;

V_dc = 100;
R = 10;
L = 0.025;
f = 50;
T = 1/f;

dt = 2*10^(-5);
t = 0:dt:5*T-dt;

Va = zeros(length(t),1);
Vb = zeros(length(t),1);
Vc = zeros(length(t),1);

Ia = zeros(length(t),1);
Ib = zeros(length(t),1);
Ic = zeros(length(t),1);

Vab = zeros(length(t),1);
Vbc = zeros(length(t),1);
Vca = zeros(length(t),1);

Van = zeros(length(t),1);
Vbn = zeros(length(t),1);
Vcn = zeros(length(t),1);

for i = 1:1:length(t)
  mymod=mod(t(i),T);
  checking=T/6;
if mymod >= 0 && mymod< checking
   Va(i) = V_dc;
   Vb(i) = 0;
   Vc(i) = V_dc;
elseif mymod>= checking  && mymod< 2*checking
   Va(i) = V_dc;
   Vb(i) = 0;
   Vc(i) = 0;
elseif mymod>= 2*checking && mymod< 3*checking
   Va(i) = V_dc;
   Vb(i) = V_dc;
   Vc(i) = 0;
elseif mymod>= 3*checking && mymod< 4*checking
   Va(i) = 0;
   Vb(i) = V_dc;
   Vc(i) = 0;
elseif mymod>= 4*checking && mymod< 5*checking
   Va(i) = 0;
   Vb(i) = V_dc;
   Vc(i) = V_dc;
elseif mymod>= 5*checking  && mymod< T
   Va(i) = 0;
   Vb(i) = 0;
   Vc(i) = V_dc;
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
     end

end

figure
% Signal names
signals = {'Va', 'Vb', 'Vc', 'Vab', 'Vbc', 'Vca', 'Van', 'Vbn', 'Vcn'};
for i = 1:length(signals)
    subplot(3, 3, i)
    plot(t, eval(signals{i}))
    ylim([-110 110])
    xlim([-0.001 0.08])
    ylabel('V(V)')
    xlabel('t(sec)')
    title(signals{i})
end
figure
% Signal names
signals = {'Ia', 'Ib', 'Ic'};
for i = 1:length(signals)
    subplot(3, 1, i)
    plot(t, eval(signals{i}))
    ylabel('I(A)')
    xlabel('t(sec)')
    title(signals{i})
end
% Calculate RMS values
Van_rms = rms(Van(1:1000));
Ia_rms = rms(Ia(1:1000));
Vbn_rms = rms(Vbn(1:1000));
Ib_rms = rms(Ib(1:1000));
Vcn_rms = rms(Vcn(1:1000));
Ic_rms = rms(Ic(1:1000));
% Calculate power output
poutput = Ia .* Van + Ib .* Vbn + Ic .* Vcn;
% Calculate apparent power
S = Van_rms * Ia_rms; % Any phase can be used
% Calculate power input
P_in = (Ia_rms .* Ia_rms) * R * ones(1000, 1); % 1000 samples per period : T/dt = 1000
% Calculate DC current
idc = sqrt(poutput / (sqrt(R^2 + (2 * pi * f * L)^2)));
figure;
subplot(2,1,1)
plot(t, idc,'g', 'LineWidth', 1.5);
ylabel('I (A)');
title('Current for Dc');
xlabel('t (sec)');
subplot(2,1,2)
plot(t, poutput,'black', 'LineWidth', 1.5);
title('Power output ');
ylabel('P (W)');
xlabel('t (sec)');
% Add a horizontal line at y = 390.21
yValue = 135;
xLimits = xlim;  % Get the x-axis limits
line(xLimits, [yValue, yValue], 'Color', 'blue', 'LineStyle', '--');

%% power factor
powerfactor = P_in(1)/S
%%All the harmonics
NFFT = 1000 ;
frequency = (1/dt)/2*linspace(0,1,NFFT/2+1);
harmonics_signals = {Van, Vbn, Vcn, Vab, Vbc, Vca};
harmonics_titles = {'Van', 'Vbn', 'Vcn', 'Vab', 'Vbc', 'Vca'};
figure;
for i = 1:numel(harmonics_signals)
    signal = harmonics_signals{i};
    title_text = [' Harmonic for ' harmonics_titles{i}];

    Vharm = fft(signal, NFFT) / 1000;

    subplot(2, 3, i);
    plot(frequency, 2*abs(Vharm(1:NFFT/2+1)));
    title(title_text);
    xlabel('F(Hz)');
end



