clear;      % clear workspace
close all;  % close all figures
clc;        % clear command window
pkg load control;   % load control package
pkg load signal;    % load signal package

%% Initialization
Vs = 300;   % input voltage
C = 10^(-5);    % capacitance
f = 10^4;   % switching frequency
R  = 3;     % resistance
dt  = 10^(-6);  % time step
needed_dimension=3; % dimension of the state space model
Lf = 5*10^(-4); % inductance
E  = 90;    % voltage across the inductor


t = 0:dt:0.05;  % simulation time
tri_wave = sawtooth(2*pi*f*t, 0.5); % generate triangular wave with amplitude 1
tri_wave = (tri_wave + 1) / 2; % shift and scale the wave to have non-negative values
plot(t, tri_wave);
hold on;
line([0,1],[0.5,0.5],'Color','r');
line([0,1],[0.8,0.8],'Color','g');
hold off;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Triangular Waveform');
xlim([0, 0.002]); % limit x-axis to 0-0.01 sec
%% 3 Cases of Duty Cycle, L
for k = 1:3
    if (k == 1)
        Dcycle = 0.5;
        L = 10^(-3);    % inductance

    elseif (k == 2)
        Dcycle = 0.5;
        L = 10^(-2);

    else
        Dcycle = 0.8;
        L = 10^(-3);
    end

% Fill Vectors for Phase 1.
% State space model of the 1st phase
A1 = [0 0 -1/L;      0 -R/Lf 1/Lf;      1/C -1/C 0];
B1 = [1/L 0;      0 -1/Lf;      0 0];
C1 = [1 0 0;       0 1 0;       0 0 1];
D1 = [0 0;         0 0;         0 0];

% Continuous-time state space model
Hcontinuous1 = ss(A1, B1, C1, D1);
 % Discretized state space model
Hdiscrete1 = c2d(Hcontinuous1, dt);

% State matrices of the discretized model
A1disc = Hdiscrete1.A;
B1disc = Hdiscrete1.B;
C1disc = Hdiscrete1.C;
D1disc = Hdiscrete1.D;

% Display the discrete A, B, C, D matrices
disp(['Discrete time matrixes for Phase 1 and  Duty cycle=:',num2str(Dcycle),'and L=',num2str(L),'.']);
disp("Discrete A matrix:");
disp(A1disc);
disp("Discrete B matrix:");
disp(B1disc);
disp("Discrete C matrix:");
disp(C1disc);
disp("Discrete D matrix:");
disp(D1disc);

% Fill Vectors for Phase 2.
% State space model of the 2nd phase
A2 = [0 0 -1/L;      0 -R/Lf 1/Lf;      1/C -1/C 0];
B2 = [0 0;          0 -1/Lf;          0 0];
C2 = [1 0 0;         0 1 0;           0 0 1];
D2 = [0 0;           0 0;             0 0];
% Continuous-time state space model
Hcontinuous2 = ss(A2, B2, C2, D2);
% Discretized state space model
Hdiscrete2 = c2d(Hcontinuous2, dt);
% State matrices of the discretized model
A2disc = Hdiscrete2.A;
B2disc = Hdiscrete2.B;
C2disc = Hdiscrete2.C;
D2disc = Hdiscrete2.D;
% Display the discrete A, B, C, D matrices
disp(['Discrete time matrixes for Phase 2 and  Duty cycle=:',num2str(Dcycle),'and L=',num2str(L),'.']);
disp("Discrete A matrix:");
disp(A2disc);
disp("Discrete B matrix:");
disp(B2disc);
disp("Discrete C matrix:");
disp(C2disc);
disp("Discrete D matrix:");
disp(D2disc);
%% Calculations
X = zeros(needed_dimension, length(t));
Y = zeros(needed_dimension, length(t));

i = 1;
while i <= length(t)
    if tri_wave(i) < Dcycle
        %complete vectors phase 1
        X(:, i+1) = A1disc*X(:, i) + B1disc*([Vs; E]);
        Y(:, i) = C1disc*X(:, i) + D1disc*([Vs; E]);
    else
        %complete vectors phase 2
        X(:, i+1) = A2disc*X(:, i) + B2disc*([Vs; E]);
        Y(:, i) = C2disc*X(:, i) + D2disc*([Vs; E]);
    end
    i = i+1;
end
% Fill plot Vectors.
Il = Y(1, :);
I_f = Y(2, :);
Vc = Y(3, :);
Ic=Y(1, :)-Y(2, :);

figure();
subplot(2,2,1);
plot(t,Il);
xlabel('time (sec)');
ylabel('I_l (A)');
title(['L Current(If) DCycle =',num2str(Dcycle),' and L = ',num2str(L)]);
xlim([0, 0.01]);

subplot(2,2,2);
plot(t,I_f);
xlabel('time (sec)');
ylabel('I_f (A)');
title(['Load Current(If) DCycle =',num2str(Dcycle),' and L = ',num2str(L)]);
xlim([0, 0.01]);


subplot(2,2,3);
plot(t,Vc);
xlabel('time (sec)');
ylabel('Vc (V)');
title(['Vc DCycle =',num2str(Dcycle),' and L = ',num2str(L)]);
xlim([0, 0.01]);


subplot(2,2,4);
plot(t,Ic);
xlabel('time (sec)');
ylabel('I_c (A)');
title(['C Current(If) DCycle =',num2str(Dcycle),' and L = ',num2str(L)]);
xlim([0, 0.01]);


end
