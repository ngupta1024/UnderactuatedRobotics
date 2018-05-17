function q1d
close all; clc; clear all;
global dt T N
N=100;
dt=0.1;
T=10;

desiredColl=q1a;
nominalColl=q1b;
xD=desiredColl(:,2:3);
uD=zeros(size(xD,1),1);
xO=nominalColl(:,2:3);
uO=nominalColl(:,1);
timestamp=linspace(0,10,100);
xFinal=xD(end,:)';%%????????
xInit=[0;0];

%set params
Qf=8*eye(2);
Q=8*eye(2);
R=1;
B=[0;1];

s2T=Qf;%2x2

% s0T=xFinal'*Qf*xFinal; %1x1
S2T=reshape(s2T,[],1); 

iterate=1;
count=0;

close all
phase=figure(1);
plot(xD(:,1),xD(:,2),'g');
hold on;
plot(xO(:,1),xO(:,2));
hold off;
xlabel('theta');
ylabel('dtheta');
title('Phase Trajectory');
xlim([-pi,pi+1]);
ylim([-2,4]);
legend('desired','nominal');
control=figure(2);
plot(timestamp,uO);
hold on;
plot(timestamp,uD);
hold off;
xlabel('t');
ylabel('u');
title('Inputs');
legend('nominal','desired');


while iterate
    count=count+1
    %use ode45 to solve Ricatti equation back in time
    [t,S]=ode45(@(t,S)ricattiEqn1(t,S,B,xO,uO,xD,uD,timestamp,Q,R), linspace(T,0,N), S2T);
    tS2=flip(t);
    S2=flipud(S);
    xBarD = xD - xO;
    s1T=-2*Qf*xBarD(end,:)';%2x1
    S1T=reshape(s1T,[],1);
    [t,S]=ode45(@(t,S)ricattiEqn2(t,S,tS2,S2,B,xO,uO,xD,uD,timestamp,Q,R), linspace(T,0,N), S1T);
    tS1=flip(t);
    S1=flipud(S);
    for i = 1:N
        s2=reshape(S2(i,:),2,2);
        s1=reshape(S1(i,:),2,1);
        K1(i,:) = (R^(-1)*B')*s2;
        K2(i)=(R^(-1)*B')*0.5*s1;
    end
    %
    [t,x] = ode45(@(t,x)dynamics(t,x,xO,uO,xD,uD,timestamp,K1,K2,tS1,tS2),timestamp,xInit);
    
    for i=1:N
        u(i) = (uD(i)) - K1(i,:)*(x(i,:)-xO(i,:))' - K2(i);
    end
    figure(control)
    hold on
    plot(t,u);
    figure(phase)
    hold on
    plot(x(:,1),x(:,2));
    
    
    drawnow;
    if norm(x(:,1)-xD(:,1))<1 && norm(x(:,2)-xD(:,2))<1.5
        fprintf('success\n');
        iterate=0;
    elseif count>=14
        fprintf('failure\n');
        iterate=0;
    else
        fprintf('thetaError=%d \n' , norm(x(:,1)-xD(:,1)));
        fprintf('thetaDError=%d \n' , norm(x(:,2)-xD(:,2)));
        xO=x;
        uO=u';
    end
end

    function xdot=dynamics(t,x,xO,uO,xD,uD,timestamp,K1,K2,tS1,tS2)
        xOInterp(1)=interp1(timestamp',xO(:,1),t);
        xOInterp(2)=interp1(timestamp',xO(:,2),t);
        xDInterp(1)=interp1(timestamp',xD(:,1),t);
        xDInterp(2)=interp1(timestamp',xD(:,2),t);
        uOInterp=interp1(timestamp',uO,t);
        uDInterp=interp1(timestamp',uD,t);
        K1Interp(1)=interp1(tS2,K1(:,1),t);
        K1Interp(2)=interp1(tS2,K1(:,2),t);
        K2Interp=interp1(tS1,K2',t);
        xBar=x-xOInterp';
        uDBar=uDInterp-uOInterp;
        u_=uDInterp-K1Interp*xBar-K2Interp;
        %     if uStar(end)>1
        %         uStar(end)=1;
        %     elseif uStar(end)<-1
        %         uStar(end)=-1;
        %     end
        xdot=[x(2); u_ - sin(x(1))-0.1*x(2)];
    end

    function Sdot=ricattiEqn1(t,S,B,xO,uO,xD,uD,timestamp,Q,R)
        xOInterp(1)=interp1(timestamp',xO(:,1),t);
        xOInterp(2)=interp1(timestamp',xO(:,2),t);
        xDInterp(1)=interp1(timestamp',xD(:,1),t);
        xDInterp(2)=interp1(timestamp',xD(:,2),t);
        uOInterp=interp1(timestamp',uO,t);
        uDInterp=interp1(timestamp',uD,t);
        A=[0 1; -cos(xOInterp(1)) -0.3];
        xbarD=xDInterp'-xOInterp';
        ubarD=uDInterp-uOInterp;
        s2_=reshape(S,2,2);
        Sdot=reshape((-(Q-s2_*B*R^(-1)*B'*s2_+s2_*A+A'*s2_)),[],1);
    end

    function Sdot=ricattiEqn2(t,S,tS2, S2,B,xO,uO,xD,uD,timestamp,Q,R)
        xOInterp(1)=interp1(timestamp',xO(:,1),t);
        xOInterp(2)=interp1(timestamp',xO(:,2),t);
        xDInterp(1)=interp1(timestamp',xD(:,1),t);
        xDInterp(2)=interp1(timestamp',xD(:,2),t);
        uOInterp=interp1(timestamp',uO,t);
        uDInterp=interp1(timestamp',uD,t);
        S2Interp(1)=interp1(tS2,S2(:,1),t);
        S2Interp(2)=interp1(tS2,S2(:,2),t);
        S2Interp(3)=interp1(tS2,S2(:,3),t);
        S2Interp(4)=interp1(tS2,S2(:,4),t);        
        A=[0 1; -cos(xOInterp(1)) -0.3];
        xbarD=xDInterp'-xOInterp';
        ubarD=uDInterp-uOInterp;
        S1_=reshape(S,2,1);
        S2=reshape(S2Interp,2,2);
        Sdot=reshape((-(-2*Q*xbarD+(A'-S2*B*R^(-1)*B')*S1_+2*S2*B*ubarD)),[],1);
    end
end