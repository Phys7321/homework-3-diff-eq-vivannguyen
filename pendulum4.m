function [kin,pot,energy,period,t,w] = pendulum4(omega0,omega,theta0,thetad0,gamma,grph) 
% Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0
    error('Must input length and initial conditions')
end

g=9.81;
%omega = sqrt(g/R); user will input omega
R = g/omega0.^2;
T= 2*pi/omega0;
% number of oscillations to graph
N = 20;


tspan = [N*T./2 N*T];
%opts = odeset('events',@events,'refine',6);
opts = odeset('refine',6);
r0 = [theta0 thetad0];

[t,w] = ode45(@proj,tspan,r0,opts,omega0,omega,gamma);
sol = [t,w];

%finding the period
ind= find(w(:,2).*circshift(w(:,2), [-1 0]) <= 0);
ind = chop(ind,4);
period= 2*mean(diff(t(ind)));

%calculating the kinetic and potential energies
kin =0.5*R^2.*sol(:,3).*sol(:,3);
pot =g*R*(1-cos(sol(:,2)));
energy = kin + pot;
%energy0 = 0.5*(R*thetad0).^2 + g*R*(1-cos(theta0));

%dE = (energy(:) - energy0)./energy0;

% Small-angle approximation solution
%delta = atan(theta0/(omega*thetad0));
%y = theta0*sin(omega*t+delta);

if grph == 1
    figure
    subplot(2,1,1)
    plot(t,w(:,1))
    title(['Position for \gamma = ' num2str(gamma)])
    xlabel('t')
    ylabel('rad')
    
    subplot(2,1,2)
    plot(t,w(:,2))
    title(['Velocity for \gamma = ' num2str(gamma)])
    xlabel('t')
    ylabel('rad/s')
    
    figure
    plot(w(:,1),w(:,2))
    title('Position vs Velocity')
    xlabel('rad')
    ylabel('rad/s')
    
if grph == 2    
    figure
    plot(w(:,1),w(:,2))
    title('Phase space')
    xlabel('theta')
    ylabel('thetad')

if grph == 3

end   
end
end
end
%-------------------------------------------
%
function rdot = proj(t,r,omega0,omega,gamma)
    rdot = [r(2); -omega0.^2*sin(r(1))-gamma.*r(2)+cos(omega.*t)];
end
