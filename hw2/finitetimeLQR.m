function u=finitetimeLQR()
close all
A=[0 1; 0 0];
B=[0;1];
x0=[0;0];
xT=[10;0];

%parameters to change for every question
Q=eye(2);
R=1;
Qf=eye(2)*100;
T=8;
initialCond= Qf;

%use ode45 to solve Ricatti equation back in time
[t_,S]=ode45(@(t,S)ricattiEqn(t,S,A,B,Q,R), [T 0], initialCond);

%plugging u back into the dynamcis and solving its ODE
%do interpolation to find S(t)
[t,x]=ode45(@(t,x)dynamics(t,x,A,B,R,linearInter(t,flipud(t_),flipud(S))),[0 T], x0);

%plot the trajectory
plot(x(:,1),x(:,2))
ax=gca;
ax.XLim=[0 10];
% ax.YLim=[-2 5];
xlabel('x');ylabel('xdot');
title('phase trajectory');

figure

for iter=1:length(t)
    S(iter,:)=linearInter(t(iter),flipud(t_),flipud(S));
    S_=reshape(S(iter,:),2,2);
    control(iter)=(-R\B'*S_)*(x(iter,:)'-[10;0]);
end
plot(control(1:end));
xlabel('time');ylabel('u');
set(gca,'XTick',[1:10:91])
set(gca,'XTickLabel',string(round(t(1:10:length(t)),4)))
title('control signal');

    function dsdt=ricattiEqn(t,S,A,B, Q, R)
        S_=reshape(S,2,2);
        dsdt_=-A'*S_-S_*A+S_*B*(R\B')*S_-Q;
        dsdt=reshape(dsdt_,4,1);
    end
    
    function xdot=dynamics(t,x,A,B,R,S)
        S_=reshape(S,2,2);
        u=(-R\B'*S_)*(x-[10;0]);
        xdot=A*(x-[10;0])+B*u;%u=-kx where u=(-R/B'*S_)[x-10;0]
    end

    function SInterp=linearInter(t,t_,S)
        if find(t_==t)
            SInterp=S(find(t_==t),:);
        else
            smaller=find(t_<t);
            greater=find(t_>t);
            
            ti=smaller(end);
            tiplus=greater(1);
            
            if abs(ti-t)<abs(tiplus-t)
                SInterp=S(ti,:);
            else
                SInterp=S(tiplus,:);
            end
        end
    end
end