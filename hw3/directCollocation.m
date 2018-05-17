function [coll]=directCollocation
%trajectory generation for a simple pendulum going 
%from 0 to pi (vertical-unstable).
close all
global dt T N
N=100;%number of collocation points
dt=0.1;%dt=T/N
T=10;%total time
xInit=[0.3;-0.3];
thresh=0.85;
time=linspace(0,T,N);
u=zeros(N,1);
guessXcoll=zeros(N,2);
initialColl=[u guessXcoll];
objective=@(initialColl) initialColl(1:N,1)'*initialColl(1:N,1); %dot(u.T,u)
A=[];
b=[];
Aeq=[];
beq=[];
ub=[repmat(thresh,N,1) inf(N,1) inf(N,1)];
lb=[repmat(-1*thresh,N,1) -inf(N,1) -inf(N,1)];
nonlinconnenq=@nonlincst;

options =optimoptions(@fmincon,'TolFun', 0.00000001,'MaxIter', 10000, ...
    'MaxFunEvals', 100000,'Display','iter', ...
    'DiffMinChange', 0.001,'Algorithm', 'sqp');
[coll,fval,ef,op]=fmincon(objective, initialColl, A,b,Aeq,beq, lb,ub,nonlinconnenq,options);

    function [c,ceq]=nonlincst(initialColl)
        xNext=initialColl(:,2:3)';
        u=initialColl(:,1);
        BoundaryFinal=[pi-xNext(1,end);0-xNext(2,end)];
        BoundaryInit=xInit-xNext(:,1);
        delta=getDefects(xNext,u');
        c = [];  %No inequality constraints
        ceq = [BoundaryInit; reshape(delta,2*N-2,1); BoundaryFinal];
    end

    %hermite-simpson method
    function delta=getDefects(xNext,u)
        xLow = xNext(:,1:(end-1));
        xHi = xNext(:,2:end);
        xdot = dynamics(xNext,u);  %no time dependence
        xdotLow = xdot(:,1:(end-1));
        xdotHi = xdot(:,2:end);
        uLow=u(:,1:end-1);
        uHi=u(:,2:end);
        xck = 0.5*(xLow + xHi)+ (dt/8).*(xdotLow-xdotHi);
        uck=0.5*(uLow+uHi);
        xdotck=[xck(2,:);-sin(xck(1,:))+uck];
        delta = xLow-xHi + (dt/6)*(xdotLow + 4*xdotck + xdotHi);
        %         delta=zeros(2,99);
    end

    function dxdt=dynamics(x,u)
        dxdt=[x(2,:);-sin(x(1,:))+u];
    end

plot(time,coll(:,1))
ylim([-1.2*thresh,1.2*thresh]);
xlabel('time');
ylabel('control input');
title('control');
figure;
plot(coll(:,2),coll(:,3))
xlabel('theta');
ylabel('dtheta');
title('Phase Trajectory');
xlim([-pi,pi+1]);
ylim([-2,2]);
end