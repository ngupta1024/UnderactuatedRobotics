function plotMat=q2a(w0)
close all
% RimlessWheel
%   rimlessWheel(w0) starts the wheel with w0 as the initial
%   rotation velocity.
%
%   rimlessWheel with no arguments starts the wheel with a random
%   initial velocity.
%
% Written by Russ Tedrake (russt@mit.edu)


m = 1; l = 1; g = 9.8; alpha = pi/8;
%gamma = 0.03;  % standing is only fixed point
gamma = 0.08;  % standing and rolling fixed points
%gamma = alpha+0.01;  % only rolling fixed point
plotMat=[];
for w0=-3:0.05:3
    if w0==0
        plotMat=[plotMat; w0 0];
    else
    x = [-sign(w0)*alpha+gamma; w0; 0];  % [\theta, \dot\theta, xfoot] 
      % I'm only keeping xfoot around so that the drawing looks smoother

    plant_dt = 1e-3;
    T = 20;
    
    for t=0:plant_dt:T
      if sign(x(2))*(x(1)-gamma) >= alpha % collision
        x = collision(x);
        plotMat=[plotMat; w0 x(2)];
        break;
      end
      x = x + plant_dt*stance_dynamics(x);  
    end
    end
end

plot(plotMat(:,1),plotMat(:,2),'.');
hold on;
plot(plotMat(:,1),plotMat(:,1));
title('Limit cycle Trajectory')
xlabel('thetaDot N');
ylabel('thetaDot N+1');
  function xp = collision(xm)

    xp = [-sign(xm(1)-gamma)*alpha + gamma; ...
      xm(2)*cos(2*alpha); ...
      xm(3) + sign(x(1)-gamma)*2*l*sin(alpha)];

  end

  function xdot = stance_dynamics(x)

    xdot = [x(2); g*sin(x(1))/l; 0];

  end
end
