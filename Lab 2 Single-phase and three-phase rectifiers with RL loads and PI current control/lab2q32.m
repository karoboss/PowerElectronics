%%Starting commands
close all;
clear;
clc;
%Starting values
V = 230;
R = 2.5;
f = 50;
%%Very helpful
cycle_period=1/f;
dt = 0.00002;
T = 10;
t = 0:dt:T-dt;
%Starting signals
Vm = 230*sqrt(2)*sqrt(3);
Vac = Vm*sin(2*pi*f*t);
Vab = Vm*sin(2*pi*f*t + pi/3);
input_length=length(Vab);
Vbc = Vm*sin(2*pi*f*t - pi/3);
Vba = -Vab;
Vca = -Vac;
Vcb = -Vbc;
%Creating the vectors frist of 0s
v_output = zeros(1,input_length);
i_output = zeros(1,input_length);

current_thyr1 = zeros(1,input_length);
current_thyr2 = zeros(1,input_length);
current_thyr3 = zeros(1,input_length);
current_thyr4 = zeros(1,input_length);
current_thyr5 = zeros(1,input_length);
current_thyr6 = zeros(1,input_length);





for L = 0.04:0.04:0.08


for a = 67

timeVecA = (a/360)*cycle_period; % seconds



            for i = 1:1:T*(1/dt)

            if mod(t(i),(cycle_period))>=timeVecA -(cycle_period)/6 && mod(t(i),(cycle_period))< timeVecA
            v_output(i) = Vcb(i);
            current_thyr1(i)=0;
            current_thyr2(i)=0;
            current_thyr3(i)=0;
            current_thyr4(i)=0;
            current_thyr5(i)=1;
            current_thyr6(i)=1;
            elseif mod(t(i),(cycle_period))>=timeVecA  && mod(t(i),(cycle_period))< timeVecA +(cycle_period)/6
            v_output(i) = Vab(i);
            current_thyr1(i)=1;
            current_thyr2(i)=0;
            current_thyr3(i)=0;
            current_thyr4(i)=0;
            current_thyr5(i)=0;
            current_thyr6(i)=1;
            elseif mod(t(i),(cycle_period))>=timeVecA +(cycle_period)/6 && mod(t(i),(cycle_period))< timeVecA +(cycle_period)/3
            v_output(i) = Vac(i);
            current_thyr1(i)=1;
            current_thyr2(i)=1;
            current_thyr3(i)=0;
            current_thyr4(i)=0;
            current_thyr5(i)=0;
            current_thyr6(i)=0;
            elseif mod(t(i),(cycle_period))>=timeVecA +(cycle_period)/3 && mod(t(i),(cycle_period))< timeVecA +(cycle_period)/2
            v_output(i) = Vbc(i);
            current_thyr1(i)=0;
            current_thyr2(i)=1;
            current_thyr3(i)=1;
            current_thyr4(i)=0;
            current_thyr5(i)=0;
            current_thyr6(i)=0;
            elseif mod(t(i),(cycle_period))>=timeVecA +(cycle_period)/2 && mod(t(i),(cycle_period))< timeVecA + (2*(cycle_period))/3
            v_output(i) = Vba(i);
            current_thyr1(i)=0;
            current_thyr2(i)=0;
            current_thyr3(i)=1;
            current_thyr4(i)=1;
            current_thyr5(i)=0;
            current_thyr6(i)=0;
            elseif mod(t(i),(cycle_period))>=timeVecA + (2*(cycle_period))/3 || mod(t(i),(cycle_period))< timeVecA - (cycle_period)/6
            v_output(i) = Vca(i);
            current_thyr1(i)=0;
            current_thyr2(i)=0;
            current_thyr3(i)=0;
            current_thyr4(i)=1;
            current_thyr5(i)=1;
            current_thyr6(i)=0;
            end

            if i > 1
            i_output(i)=L/(R*dt+L)*i_output(i-1) + dt/(R*dt+L)*v_output(i);
            end

            if  i_output(i) <= 0
                i_output(i)=0;
            end

            current_thyr1(i) = current_thyr1(i)*i_output(i);
            current_thyr2(i) = current_thyr2(i)*i_output(i);
            current_thyr3(i) = current_thyr3(i)*i_output(i);
            current_thyr4(i) = current_thyr4(i)*i_output(i);
            current_thyr5(i) = current_thyr5(i)*i_output(i);
            current_thyr6(i) = current_thyr6(i)*i_output(i);



            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING Vin Vout Iout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
grid on ;
subplot(2,1,1)
plot(t,Vac);
grid on;
hold on;
plot(t,Vab);
hold on;
plot(t,Vbc);
hold on;
plot(t,Vba);
hold on;
plot(t,Vca);
hold on;
plot(t,Vcb);
hold on;
plot(t,v_output,'green','LineWidth', 3);
xlim([3 3.04])
ylabel('Vinput(V)')
xlabel('time(sec)')
title(['Voutput for  a = ',num2str(a),' and L = ',num2str(L),','])
subplot(2,1,2)
plot(t,i_output,'red','LineWidth', 3)
grid on
ylabel('Ioutput(A)')
xlabel('t(sec)')
xlim([3 3.02])
title(['Ioutput for  a = ',num2str(a),' and L = ',num2str(L),','])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING THYRISTOR CURRENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
grid on
hold on;
plot(t,i_output,'green','LineWidth', 4);
hold on;
grid on;
plot(t,current_thyr1);
hold on;
plot(t,current_thyr2);
hold on;
plot(t,current_thyr3);
hold on;
plot(t,current_thyr4);
hold on;
plot(t,current_thyr5);
hold on;
plot(t,current_thyr6);
hold off;
legend('Io','CurrentThyr1','CurrentThyr2','CurrentThyr3','CurrentThyr4','CurrentThyr5','CurrentThyr6');
ylabel('I(A)');
xlabel('t(sec)');
xlim([3 3.04]);
ylim([0 300]);
title('Thyristor currents plots')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING ALL TOGETHER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
grid on
hold on;
plot(t,v_output,'red','LineWidth', 4);
plot(t,Vac);
grid on;
hold on;
plot(t,Vab);
hold on;
plot(t,Vbc);
hold on;
plot(t,Vba);
hold on;
plot(t,Vca);
hold on;
plot(t,Vcb);
plot(t,i_output,'green','LineWidth', 4);
hold on;
grid on;
plot(t,current_thyr1);
hold on;
plot(t,current_thyr2);
hold on;
plot(t,current_thyr3);
hold on;
plot(t,current_thyr4);
hold on;
plot(t,current_thyr5);
hold on;
plot(t,current_thyr6);
hold off;
legend('Vo','Io','CurrentThyr1','CurrentThyr2','CurrentThyr3','CurrentThyr4','CurrentThyr5','CurrentThyr6');
ylabel('Values');
xlabel('t(sec)');
xlim([3 3.04]);
title('Combined all together plot')



end
end
