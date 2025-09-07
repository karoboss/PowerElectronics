%%Starting commands
close all;
clear;
clc;
%Starting values
V = 230;
R = 2.5;
f = 50;
dt = 0.00002;
T = 10;
t = 0:dt:T;
cycle_period=1/f;
%%Time until 3 seconds
t1to3 = 1:round(2*length(t)/5);
%The part we want t=5 s from 4 to 9
t4to9 = 1+round(2*length(t)/5):round(9*length(t)/10);
%last 1 seconds
t9to10 = round(9*length(t)/10)+1:length(t);

Vm = V * sqrt(2);
Vm = Vm * sqrt(3);
v_down = 200 * sqrt(2);
v_down = v_down * sqrt(3);

%3 First seconds with 230V
Vab(t1to3) = Vm * sin((2*pi*f).*t(t1to3));
Vac(t1to3) = Vm * sin((2*pi*f).*t(t1to3) - pi/3);
Vbc(t1to3) = Vm * sin((2*pi*f).*t(t1to3) - 2*pi/3);

%Drop for 5 secs volts to 200
Vab(t4to9) = v_down*sin((2*pi*f).*t(t4to9));
Vac(t4to9) = v_down*sin((2*pi*f).*t(t4to9) - pi/3);
Vbc(t4to9) = v_down*sin((2*pi*f).*t(t4to9) - 2*pi/3);

%1 last seconds 230V
Vab(t9to10) = Vm*sin((2*pi*f).*t(t9to10));
Vac(t9to10) = Vm*sin((2*pi*f).*t(t9to10) - pi/3);
Vbc(t9to10) = Vm*sin((2*pi*f).*t(t9to10) - 2*pi/3);
input_length=length(Vab);
Vba = -Vab;
Vca = -Vac;
Vcb = -Vbc;
%From the excersice
i_mean = 75;
t_sampling = 0.02;
test=t(1001:round(t_sampling/dt-2):end);

%Creating the vectors frist of 0s
v_output = zeros(1,input_length);
i_output = zeros(1,input_length);

current_thyr1 = zeros(1,input_length);
current_thyr2 = zeros(1,input_length);
current_thyr3 = zeros(1,input_length);
current_thyr4 = zeros(1,input_length);
current_thyr5 = zeros(1,input_length);
current_thyr6 = zeros(1,input_length);

i_start = zeros(1,10/t_sampling);
a_start = zeros(1,10/t_sampling);

temp = 2;
a = 0;
L = 0.04;

for i = 1:input_length
%Phase in angle
    phase = mod(t(i),(cycle_period))*360/(cycle_period);
%a adjustment
    if t(i)>=t_sampling && mod(t(i),t_sampling) == 0

        t1 = round((t(i)-t_sampling)/dt)+1;
        t2 = round(t(i)/dt)+1;

        i_start(temp) = mean(i_output(t1:t2));

        % settign my Po
        e0 = i_start(temp-1)-i_mean;
        e1 = i_start(temp)-i_mean;
        P = 0.1;
        I = 7.9;
        a = a + P*(e1-e0) + I*e1*t_sampling;
        if a<0
        	a = 0;
        end

        if a>90
            a = 90;
        end

        a_start(temp) = a;
        temp = temp +1;

    end
    if a <= 60
        if  phase >= a && phase < a + 60
            v_output(i) = Vcb(i);
            current_thyr5(i) = 1;
            current_thyr6(i) = 1;
        elseif phase >= a + 60 && phase < a + 120
            v_output(i) = Vab(i);
            current_thyr1(i) = 1;
            current_thyr6(i) = 1;
        elseif phase >= a + 120 && phase < a + 180
            v_output(i) = Vac(i);
            current_thyr1(i) = 1;
            current_thyr2(i) = 1;
        elseif phase >= a + 180 && phase < a + 240
            v_output(i) = Vbc(i);
            current_thyr3(i) = 1;
            current_thyr2(i) = 1;
        elseif phase >= a + 240 && phase < a + 300
            v_output(i) = Vba(i);
            current_thyr3(i) = 1;
            current_thyr4(i) = 1;
        elseif phase >= a + 300 || phase < a
            v_output(i) = Vca(i);
            current_thyr5(i) = 1;
            current_thyr4(i) = 1;
        end
    else %for a>60
        if phase >= a - 60 && phase < a
            v_output(i) = Vca(i);
        	current_thyr5(i) = 1;
            current_thyr4(i) = 1;
        elseif phase >= a && phase < a + 60
            v_output(i) = Vcb(i);
            current_thyr5(i) = 1;
            current_thyr6(i) = 1;
        elseif phase >= a + 60 && phase < a + 120
            v_output(i) = Vab(i);
            current_thyr1(i) = 1;
            current_thyr6(i) = 1;
        elseif phase >= a + 120 && phase < a + 180
            v_output(i) = Vac(i);
            current_thyr1(i) = 1;
            current_thyr2(i) = 1;
        elseif phase >= a + 180 && phase < a + 240
            v_output(i) = Vbc(i);
            current_thyr3(i) = 1;
            current_thyr2(i) = 1;
        elseif  phase >= a + 240 || phase < a - 60
            v_output(i) = Vba(i);
            current_thyr3(i) = 1;
            current_thyr4(i) = 1;
        end
    end


    if i>1
        i_output(i)=L/(R*dt+L)*i_output(i-1) + dt/(R*dt+L)*v_output(i);
    end

    if  i_output(i)<=0
        v_output(i) = 0;
        i_output(i) = 0;
    end

    current_thyr1(i) = current_thyr1(i)*i_output(i);
    current_thyr2(i) = current_thyr2(i)*i_output(i);
    current_thyr3(i) = current_thyr3(i)*i_output(i);
    current_thyr4(i) = current_thyr4(i)*i_output(i);
    current_thyr5(i) = current_thyr5(i)*i_output(i);
    current_thyr6(i) = current_thyr6(i)*i_output(i);


end

figure();
hold on

subplot(3,1,1)
plot(t,v_output,'green')
title('Voutput')
ylabel('Voutput(V)')
xlabel('time(sec)')
subplot(3,1,2)
plot(t,i_output,'magenta')
title('Ioutput')
ylabel('Ioutput(A)')
xlabel('time(sec)')
subplot(3,1,3)
plot(test,a_start)
title('a')
ylim([0 100])
ylabel('degrees')
xlabel('time(sec)')
