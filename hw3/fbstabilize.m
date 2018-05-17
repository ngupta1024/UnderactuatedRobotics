function fbstabilize
close all; clear all; clc;
global T dt N uActual
T=10;
N=100;
dt=0.1;
boundval=1;
optimalColl=directCollocation;
xOptimal=optimalColl(:,2:3);
uOptimal=optimalColl(:,1);
timestamp=linspace(0,10,100);
%set params
Qf=eye(2);
xInit=[0.3;-0.3];
Q=eye(2);
R=1;
B=[0;1];
%use ode45 to solve Ricatti equation back in time
[t,S]=ode45(@(t,S)ricattiEqn(t,S,Q,R,B,xOptimal,timestamp), linspace(T,0,N), Qf);
for i = 1:N   
    Sreshaped(:,:,i) = reshape(S(length(t)-i+1,:),2,2);
    K(i,:) = (R\B')*Sreshaped(:,:,i);
end

[t,x] = ode45(@(t,x)dynamics(t,x,xOptimal,uOptimal, timestamp,K,B),timestamp,xInit);
    
    function xdot=dynamics(t,x,xOptimal,uOptimal,timestamp,K,B)
        A=pfl(t,xOptimal,timestamp);
        KInterp=[linInterp(t,K(:,1)',timestamp) linInterp(t,K(:,2)',timestamp)];
        uExpected=linInterp(t,uOptimal',timestamp);
        xExpected=[linInterp(t,xOptimal(:,1)',timestamp);...
            linInterp(t,xOptimal(:,2)',timestamp)];
        xbar=x-xExpected;
        uStar=uExpected-KInterp*(xbar);
        if (uStar>boundval)
            uStar = boundval;
        elseif (uStar<-1*boundval)
            uStar = -1*boundval;
        end
        uActual(end+1)= uStar;
        xdot=[x(2); uStar - sin(x(1))];%u=-kx where u=(-R/B'*S_)[x-10;0]
    end

    function u = linInterp(t,u,timestamp)    
        timestamp = [timestamp 10.1];
        u = [u u(end)];
        %Find nearest points
        index = find(timestamp <= t,1,'last');
        t1 = timestamp(index);
        t2 = timestamp(index+1);
        u1 = u(index);
        u2 = u(index+1);
        %Linear Interpolation
        u = u1 + ((t-t1)*((u2-u1)/(t2-t1)));
    end

    function dsdt=ricattiEqn(t,S,Q,R,B,xOptimal,timestamp)
        A=pfl(t,xOptimal,timestamp);
        S_=reshape(S,2,2);
        dsdt_=-A'*S_-S_*A+S_*B*(R\B')*S_-Q;
        dsdt=reshape(dsdt_,4,1);
    end

    function A=pfl(t,xOpt,timestamp)
        xOpt=[xOpt; xOpt(end,:)];
        idx = find(timestamp <= t,1,'last');
        tsmod=[timestamp 10.1];
        t1 = tsmod(idx);
        t2 = tsmod(idx+1);
        x1 = xOpt(idx,1);
        x2 = xOpt(idx+1,1);
        x = x1 + ((t-t1)*((x2-x1)/(t2-t1)));
        A=[0 1; cos(x) 0];
    end
figure;
plot(x(:,1),x(:,2),'r');
hold on;
plot(xOptimal(:,1),xOptimal(:,2));
hold off;
xlabel('theta');
ylabel('dtheta');
title('Phase Trajectory');
legend('actual','optimal');
figure;
plot(linspace(0,10,length(uActual)),smooth(uActual,15));
hold on;
plot(timestamp,uOptimal);
hold off;
xlabel('t');
ylabel('u');
title('Inputs');
legend('actual','optimal');
end