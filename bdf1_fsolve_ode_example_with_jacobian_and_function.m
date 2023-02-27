%bdf1+fsolve solver for stiff equations with input Jacobian matrix and function handle:

%memory cleanup:
clear all;
close all;
clc;
tic;

%parameters:
t0=0;
tfinal=2;
h=0.1; %step size
t=t0:h:tfinal; %time vector
n=length(t); %time points
y=zeros(1,n); %solution vector
y0 = 0; % initial condition
y(1)=y0; %initiating solution vector
options = optimset('Display','off','Jacobian','on'); %fsolve optimization

% Use backward Euler method and fsolve to solve the ODE
for i=2:length(t)
    y0 = y(i-1); %previous step solution
    y(i) = fsolve( @(y) BEJ(y, y0, t(i-1),h,@ode,@ode_j), y0, options);
end

% Analytical solution plot
t1=0:0.01:2; %time vector
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
toc;

%the ODE to solve:
function f=ode(t,y)
   f= 50*(cos(t)-y);
end

%Jacobian of the ODE:
function j=ode_j()
   j= -50;
end

function [F,J] = BEJ(y, y0, t, h,ode,ode_j)
       F=y-y0-h*ode(t,y); %function input for fsolve in terms of bdf1
       J=1-h*ode_j(); %function jacobian
end
