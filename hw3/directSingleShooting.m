function directSingleShooting
clear all; close all; clc;
%%Before submitting make plots for every iteration
%%simple pendulum
global xInitial
global T
global dt
thresh=0.85;

dt=0.025;
T=200*dt;
time=0:dt:T;
u0=zeros(length(time),1);
xInitial=[0;0];
%xNext=simulateTraj(xInitial,u0,T,dt);
objective=@(u) u'*u;
A=[];
b=[];
Aeq=[];
beq=[];
ub=repmat(thresh,201,1);
lb=repmat(-1*thresh,201,1);
% [C,Ceq]=nonlincst(u);
nonlinconnenq=@nonlincst;
options =optimoptions(@fmincon,'TolFun', 0.00000001,'MaxIter', 10000, ...
    'MaxFunEvals', 100000,'Display','iter', ...
    'DiffMinChange', 0.001,'Algorithm', 'sqp');
[u,fval,ef,op]=fmincon(objective, u0, A,b,Aeq,beq, lb,ub,nonlinconnenq,options);
xNext=simulateTraj(xInitial,u,T,dt);
u(end)=u(end-1);

    function xNext=simulateTraj(xInitial,u,T,dt)
        xNext(:,1)=xInitial;
        for iter=2:T/dt+1
            x=xNext(:,iter-1);
            dx1dt=x(2);
            dx2dt=-sin(x(1))+u(iter-1);
            xNext(:,iter)=[x(1)+ dt*dx1dt; x(2)+dt*dx2dt];
        end
    end

%non linear constraint function
    function [C,Ceq]=nonlincst(u)
        xNext=simulateTraj(xInitial,u,T,dt);
        Ceq=[pi-xNext(1,end);0-xNext(2, end)];
        C=[];
    end

plot(xNext(1,:),xNext(2,:));
xlabel('theta');
ylabel('dtheta');
title('Phase Trajectory');
xlim([-pi,pi+1]);
ylim([-2,2]);
figure;
plot(time,u);
xlim([0,5]);
ylim([-1.2*thresh,1.2*thresh]);
xlabel('time');
ylabel('control input');
title('control');
end