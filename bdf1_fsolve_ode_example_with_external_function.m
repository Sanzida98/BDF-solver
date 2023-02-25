%bdf1+fsolve solver for stiff equations with function input:

%memory cleanup:
clc;
clear all;
clearvars;

%parameters:
t0=0;
tfinal=2;
h=0.1; %step size
n=ceil((tfinal-t0)/h);
t=zeros(1,n); %time vector
t(1)=t0;
y=zeros(1,n); %solution vector
y0 = 0; % initial condition
y(1)=y0; %initiating solution vector

% Function definition
F = @(y0,y,t) y-y0-h*ode_fun(t,y);  %function in terms of bdf1
options = optimset('Display','off');

% Use backward Euler method and fsolve to solve the ODE
for i=1:n
    t(i+1)=t(i)+h;
    y0 = y(i); %previous step solution
    y(i+1) = fsolve(@(y) F(y0, y, t(i)), y0,options);
end

% Analytical solution plot
t1=0:0.01:2; %time vector
y_exact= (50/2501)*(sin(t1)+50*cos(t1)-50*exp(-50*t1)); %analytical solution plot

% Numerical Solution plot
figure();
plot(t,y,'-r','LineWidth',3);
grid on;
hold on;
plot(t1, y_exact,'-b','LineWidth',3);
legend('Numerical Solution','Analytical Solution');
xlabel('time(t)');
ylabel('solution(y)');