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
%Time vector
t = 0:dt:T;
%Starting signals
v_input = V * sqrt(2)* sin((2*pi*f).*t);
input_length=length(v_input);
%Creating the vectors frist of 0s
v_output = zeros(1,input_length);
i_output= zeros(1,input_length);
current_thyr1_2 = zeros(1,input_length);
current_thyr3_4 = zeros(1,input_length);
%%Create the for loop create the Voutput and the Ioutput based on the thyristors
for a = 0:90:90;
%%This is the a in the time
timeVecA = (a/360)*cycle_period; % seconds

    for L= 0.04:0.04:0.08;
       for i = 1:input_length
          %Here we could use the bottom type but i thought it way to easier to eliminate the constante values
          %if mod((2*pi*f*t(i)),(2*pi)) >= timeVecA && mod((2*pi*f*t(i)),(2*pi)) <= timeVecA +2*pi/2
          %cycle_period=1/f 2*pi/2*pi*f


            %For the first condition we have positive supply voltage and 1and2 thyristors working
            if mod(t(i),cycle_period)>=timeVecA && mod(t(i),cycle_period)<=timeVecA + cycle_period/2
                v_output(i) = v_input(i);
                current_thyr1_2(i) = 1;
                current_thyr3_4(i) = 0;
                %Else we have negative sypply and 34 thyristors working
            else
                v_output(i) = -v_input(i);
                current_thyr1_2(i) = 0;
                current_thyr3_4(i) = 1;
            end
          %%Type to calculate the ioutput
            if i>1
                i_output(i)=L/(R*dt+L)*i_output(i-1) + dt/(R*dt+L)*v_output(i);

            end
            %Case we have negative current
            if i_output(i)<=0
                v_output(i) = 0;
                i_output(i) = 0;
            end
            %final calculations of the thyristors currents
            current_thyr1_2(i) = current_thyr1_2(i)*i_output(i);
            current_thyr3_4(i) = current_thyr3_4(i)*i_output(i);

        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING Vin Vout Iout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure();
        grid on ;
        subplot(3,1,1)
        plot(t,v_input,'-.','LineWidth', 3)
        xlim([3 3.04])
        ylabel('Vinput(V)')
        xlabel('time(sec)')
        title(['Vinput for  a = ',num2str(a),' and L = ',num2str(L),','])
        subplot(3,1,2)
        plot(t,v_output,'blue','LineWidth', 3)
        grid on
        xlim([3 3.04])
        ylabel('Voutput(V)')
        xlabel('time(sec)')
        title(['Voutput for  a = ',num2str(a),' and L = ',num2str(L),','])
        subplot(3,1,3)
        plot(t,i_output,'red','LineWidth', 3)
        grid on
        ylabel('Ioutput(a)')
        xlabel('t(sec)')
        xlim([3 3.04])
        title(['Ioutput for  a = ',num2str(a),' and L = ',num2str(L),','])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING THYRISTOR CURRENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure();
        grid on
        hold on;
        plot(t,current_thyr1_2,'-.','LineWidth', 3)
        plot(t,current_thyr3_4,'LineWidth', 3)
        xlim([3 3.04])
        legend('CurrentThyr12','CurrentThyr34')
        title('Thyristor currents plots')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING ALL TOGETHER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure();
        hold on
        grid on
        plot(t,v_input,'-.','LineWidth', 3)
        plot(t,v_output,'green','LineWidth', 3)
        plot(t,i_output,'red','LineWidth', 3)
        plot(t,current_thyr1_2,'magenta','LineWidth', 3)
        plot(t,current_thyr3_4,'black','LineWidth', 3)
        xlim([3 3.04])
        xlabel('time(sec)')
        ylabel('Values')
        legend('Vs','Voutput','Ioutput','CurrentThyr12','CurrentThyr34')
        title('Combined all together plot')

    end
end




