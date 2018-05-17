function q2b

x = [-0.3127; 4.5]; %initial guess
threshold = 0.001;
%F=x-P(x); xstar is the root
F = [inf; inf];
m = 0.1;
while (norm(F) > threshold)
    F = x - P(x);
    Fderiv = eye(2) - dpdx(x)+m*eye(2); %to prevent it from being singular
    x = x - Fderiv\F;
    fprintf('x value:%d, %d\n',x(1),x(2));
end
fixedPoint = x;

%Jacobian around the fixed point
J = dpdx(fixedPoint);
fprintf('Eigenvalues of Jacobian: %d\n ', eig(J));

    function xNext = P(x)
        %States
        theta = x(1);
        dtheta = x(2);
        
        %Parameters
        m = 1; l = 1; g = 9.8; alpha = pi/8; gamma = 0.08;
        
        w1 = sqrt((2*g/l)*(1 - cos(gamma-alpha)));
        w2 = -sqrt((2*g/l)*(1 - cos(alpha-gamma)));
        
        if (dtheta>w1)
            dthetaNext = cos(2*alpha)*sqrt((dtheta)^2 + 4*(g/l)*sin(alpha)*sin(gamma));
        elseif (dtheta<=w1) && (dtheta>=w2)
            dthetaNext = -dtheta*cos(2*alpha);
        elseif (dtheta<w2)
            dthetaNext = -cos(2*alpha)*sqrt((dtheta)^2 + 4*(g/l)*sin(alpha)*sin(gamma));
        else
            dthetaNext=dtheta;
        end
        
        xNext = [theta; dthetaNext];
    end


    function J = dpdx(x)
        theta = x(1);
        dtheta = x(2);
        
        %Parameters
        m = 1; l = 1; g = 9.8; alpha = pi/8; gamma = 0.08;
        
        w1 = sqrt((2*g/l)*(1 - cos(gamma-alpha)));
        w2 = -sqrt((2*g/l)*(1 - cos(alpha-gamma)));
        
        %found analytically
        if (dtheta>w1)
            dxNextdx = (dtheta*cos(2*alpha))/(dtheta^2 + (4*g*sin(alpha)*sin(gamma))/l)^(1/2);
        elseif (dtheta<=w1) && (dtheta>=w2)
            dxNextdx = -cos(2*alpha);
        elseif (dtheta<w2)
            dxNextdx = -(dtheta*cos(2*alpha))/(dtheta^2 + (4*g*sin(alpha)*sin(gamma))/l)^(1/2);
        else
            dxNextdx=randn;
        end
        
        J = [1 0;...
            0 dxNextdx];
    end
end

