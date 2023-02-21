%BDF+Newton's method solver for Stiff ODEs:
clc;
close all;
clearvars;
clear all;
set(0,'defaultlinelinewidth',3);
set(0,'defaultTextFontSize',18);
set(0,'DefaultAxesFontSize',24);
%parameters:
t0 = 0; %initial time
h=0.1; %step size
tfinal = 2; %final time
n = ceil((tfinal-t0)/h); %number of time points
N = 20; %Iteration points for Newton's method
t=zeros(1,n); %time vector in the range t0<=h<=tfinal
y=zeros(1,n); %ODE solution vector
t(1)  =  t0; %starting time
y0 = 0; %initial condition
y(1) =  y0; %initiating with initial condition

%function specification:
F=@(t,y,y0) y-y0-h*ode_fun(t,y); %Backward Euler function: y(i)-y(i-1)-h*f(t(i),y(i))=0
dF=@(t,y) (1-h*ode_fun_jac(t,y)); %dF/dy

%loop for backward Euler:
for i=1:n
    y0 = y(i); %last step solution update
    t(i+1)=t(i)+h;

%loop for Newton's method:
    for j=1:N
        x = y(i)-(F(t(i),y(i),y0)/dF(t(i),y(i)));
    end
    y(i+1)=x;
end
% Analytical solution plot
t1=0:0.01:2; %time vector
y_exact= (50/2501)*(sin(t1)+50*cos(t1)-50*exp(-50*t1)); %analytical solution plot
%solution plot:
figure()
plot(t,y);
hold on;
grid on;
ylabel('y(t)');
xlabel('t');
plot(t1, y_exact);
legend('Numerical Solution','Analytical Solution');
hold off;
