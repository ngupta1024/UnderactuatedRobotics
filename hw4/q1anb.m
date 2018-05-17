function q1anb
close all; clc; clear all;
global dt T N
N=100;
dt=0.1;
T=10;

desiredColl=q1a;
nominalColl=q1b;
xD=desiredColl(:,2:3);
uD=desiredColl(:,1);
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

end