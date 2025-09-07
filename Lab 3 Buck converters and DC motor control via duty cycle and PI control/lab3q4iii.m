clear all;
close all;
clc;

Vs_A = 300;
L = 10^(-3);
C = 10^(-5);
f = 10^4;
w_ref = 80;

Ra = 3;
La = 0.5*10^(-3);
J = 5*10^(-3);
kt = 0.3;
ke = 0.3;

dt = 10^(-6);

t = 0:dt:2;
tri_wave = sawtooth(2*pi*f*t, 0.5); % generate triangular wave with amplitude 1
tri_wave = (tri_wave + 1) / 2; % shift and scale the wave to have non-negative values

Il = zeros(length(t), 1);
If = zeros(length(t), 1);
Vc = zeros(length(t), 1);
Icoil= zeros(length(t), 1);
W = zeros(length(t), 1);
duty = zeros(length(t), 1);
error = zeros(length(t), 1);
X = zeros(length(t), 1);
Te = zeros(length(t), 1);
xprev=zeros(4, 1);
%Phase1
A1 = [0 0 -1/L 0;
      0 -Ra/La 1/La -ke/L;
      1/C -1/C 0 0;
      0 kt/J 0 0];

B1 = [1/L 0;
      0 0;
      0 0;
      0 -1/J];

C1 = [1 0 0 0;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];

D1 = [0 0;
      0 0;
      0 0;
      0 0];


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
disp(['Discrete time matrixes for Phase 1']);
disp("Discrete A matrix:");
disp(A1disc);
disp("Discrete B matrix:");
disp(B1disc);
disp("Discrete C matrix:");
disp(C1disc);
disp("Discrete D matrix:");
disp(D1disc);

%Phase2
A2 = [0 0 -1/L 0;
      0 -Ra/La 1/La -ke/L;
      1/C -1/C 0 0;
      0 kt/J 0 0];

B2 = [0 0;
      0 0;
      0 0;
      0 -1/J];

C2 = [1 0 0 0;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];

D2 = [0 0;
      0 0;
      0 0;
      0 0];

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
disp(['Discrete time matrixes for Phase 2 ']);
disp("Discrete A matrix:");
disp(A2disc);
disp("Discrete B matrix:");
disp(B2disc);
disp("Discrete C matrix:");
disp(C2disc);
disp("Discrete D matrix:");
disp(D2disc);
%PI Rithmisi

P=1.2;
I=13.6;

A=0;
B=1;
C=I;
D=P;

sys_PI=ss(A, B, C, D);
sys_PId=c2d(sys_PI, dt);

Ad = sys_PId.A;
Bd = sys_PId.B;
Cd = sys_PId.C;
Dd = sys_PId.D;

disp(['Discrete time matrixes for Controller ']);
disp("Discrete A matrix:");
disp(Ad);
disp("Discrete B matrix:");
disp(Bd);
disp("Discrete C matrix:");
disp(Cd);
disp("Discrete D matrix:");
disp(Dd);






    for i=1:length(t)-1 % change loop range to avoid index out of bounds
        if(t(i)>1)
      TL=25;
  else
      TL=20;
  end


       error(i) = w_ref - W(i); % update error for current step

X(i+1) = Ad * X(i) + Bd * error(i); % update controller state for current step
duty(i+1) = Cd * X(i) + Dd * error(i); % compute duty cycle for current step


if  tri_wave(i+1) < duty(i+1)/100 % use i+1 instead of i to get the current value instead of the previous
    % Vs = 300
    A_phase = A1disc;
    B_phase = B1disc;
    C_phase = C1disc;
    D_phase = D1disc;
    u = [Vs_A; TL];
else
    % Vs = 0
    A_phase = A2disc;
    B_phase = B2disc;
    C_phase = C2disc;
    D_phase = D2disc;
    u = [0; TL];
end
%%OI EKSISOSEIS MAS
x = A_phase * xprev + B_phase * u;
y = C_phase * x + D_phase * u;
%%PRINT OUR ELEMENTS
Il(i+1) = y(1);
If(i+1) = y(2);
Icoil(i+1)=y(1)-y(2);
Vc(i+1) = y(3);
W(i+1) = y(4);
Te(i+1) = kt * If(i+1);

xprev = x; % update previous state for next iteration
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%

    figure(1);
    plot(t, Il,'m');
    title(['Current coil for TL=',num2str(TL),' and 20']);
    ylabel('I');
    xlabel('time');

    figure(2);
    plot(t, If,'m');
    title(['Current load for TL=',num2str(TL),' and 20']);
    ylabel('I');
    xlabel('time');


    figure(3);
    plot(t, Icoil,'m');
    title(['Capacitor current for TL=',num2str(TL),' and 20']);
     ylabel('I');
    xlabel('time');


    figure(4);
    plot(t, Vc,'m');
    title(['Vc for TL=',num2str(TL),' and 20']);
    ylabel('V(V)');
    xlabel('time');

    figure(5);
    plot(t, W,'m');
    title(['Omega for TL=',num2str(TL),' and 20']);
    ylabel('Omega(rad/s)');
    xlabel('time');

    figure(6);
    plot(t, Te,'m');
    title(['Torque for TL=',num2str(TL),' and 20']);
    ylabel('Torque(Nm)');
    xlabel('time');



