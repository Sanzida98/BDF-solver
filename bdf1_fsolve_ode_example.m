%bdf1+fsolve solver for stiff equations:

%memory cleanup:
clear all;
close all;
clc;
set(0,'defaultlinelinewidth',3);
set(0,'defaultTextFontSize',18);
set(0,'DefaultAxesFontSize',24);
%parameters:
t0=0;
tfinal=200;
h=20; %step size
t=t0:h:tfinal; %time vector
n=length(t); %time points
y=zeros(1,n); %solution vector
y0 = 0; % initial condition
y(1)=y0; %initiating solution vector

% Function definition
F = @(y0, y, t) y-y0-h*ode_fun(t,y);  %function in terms of bdf1

% Use backward Euler method and fsolve to solve the ODE
for i=2:length(t)
    y0 = y(i-1); %previous step solution
    y(i) = fsolve(@(y) F(y0, y, t(i-1)), y0);
end

% Analytical solution plot
t1=0:20:200; %time vector
y_exact= (50/2501)*(sin(t1)+50*cos(t1)-50*exp(-50*t1)); %analytical solution plot

% Solution plot
figure();
plot(t,y,'-r','LineWidth',3);
grid on;
hold on;
plot(t1, y_exact,'-b','LineWidth',3);
legend('Numerical Solution','Analytical Solution');
xlabel('time(t)');
ylabel('solution(y)');
