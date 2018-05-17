function xTrajectory=infLQRCartPole(x)
close all

if nargin==0
    x=[0 2.4 0 0]';
end

x(2)=wrapTo2Pi(x(2));
plt=true;
T = 10;
plant_dt = 0.01;

% Linearized dynamics
A=[0 0 1 0; 0 0 0 1; 0 1 0 0; 0 2 0 0];
B=[0 0 1 1]';

% LQR
Q = 10*eye(4,4);
R = 1;

[K,S] = lqr(A,B,Q,R);

xd = [0, pi, 0, 0]';
Ed = 1;


xTrajectory=[0,x'];
for t=plant_dt:plant_dt:T
    xdiff=(x-xd);
    xdiff(2) = mod(xdiff(2)+pi, 2*pi)-pi;
    u = -K*xdiff; 
    xdot = dynamics(x,u);
    x = x + plant_dt*xdot;
    xTrajectory(end+1,:)=[t,x'];
end

if plt==true
    plot(xTrajectory(:,2),xTrajectory(:,4)); 
    xlim([-2.5,2.5]);
    ylim([-3.5,3.5]);
    title('Trajectory(x) for the border-line case')
    figure;
    plot(xTrajectory(:,3),xTrajectory(:,5)); 
    xlim([0,2*pi]);
    ylim([-2.5,2.5]);
    title('Trajectory(theta) for the border-line case')
end
    function xdot = dynamics(x,u)
        s = sin(x(2)); c = cos(x(2));
        xddot = [u + s*x(4)^2 + s*c]/[1+s^2];
        tddot = [-u*c - x(4)^2*c*s - 2*s]/[1+s^2];
        xdot = [x(3:4); xddot; tddot];
    end
end