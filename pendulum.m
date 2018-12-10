function [period,sol] = pendulum(R,theta0,thetad0,grph) 
% Finds the period of a nonlinear pendulum given the length of the pendulum
% arm and initial conditions. All angles in radians.

%Setting initial conditions
if nargin==0
    error('Must input length and initial conditions')
end
if nargin==1 %if only 1 input R, then set theta to pi/2 and thetad0 to 0 and use grph 0
   theta0 = pi/2;
   thetad0=0;
   grph=0;
end
if nargin==2 %if only 2 input R and theta0, then set thetad0 to 0 and use grph 1
    thetad0 = 0;
    grph=1;
end
if nargin==3 %if only 3 input R and theta0 and thetad0, then use grph 1
    grph=1;
end

g=9.81;
omega = sqrt(g/R);
T= 2*pi/omega;
% number of oscillations to graph
N = 10;


tspan = [0 N*T];
%opts = odeset('events',@events,'refine',6); %Here for future event finder
opts = odeset('refine',6);
r0 = [theta0 thetad0];

[t,w] = ode45(@proj,tspan,r0,opts,g,R);
sol = [t,w];

%finding the period
ind= find(w(:,2).*circshift(w(:,2), [-1 0]) <= 0);
ind = chop(ind,4);
period= 2*mean(diff(t(ind)));

%calculating the potential and kinetic energies
pot = g*R*(1-cos(w(ind(1):ind(3),1)));
kin =  0.5*(R*w(ind(1):ind(3),2)).^2;
energy = pot + kin;
energy0 = 0.5*(R*thetad0).^2 + g*R*(1-cos(theta0));

dE = (energy(:) - energy0)./energy0;

% Small-angle approximation solution
delta = atan(theta0/(omega*thetad0));
y = theta0*sin(omega*t+delta);

if grph % Plot Solutions of exact and small angle
    figure
    subplot(2,1,1)
    plot(t,w(:,1),'k-',t,y,'b--')
    legend('Exact','Small Angle')
    title('Exact vs Approximate Solutions')
    xlabel('t')
    ylabel('\phi')
    subplot(2,1,2)
    plot(t,w(:,1)-y,'k-')
    title('Difference between Exact and Approximate')
    xlabel('t')
    ylabel( '\Delta\phi')
    
    figure
    plot(t(ind(1):ind(3)),dE,'m:')
    title('Total Energy during one cycle')
    xlabel('t')
    ylabel('dn')
    
    figure
    subplot(2,1,1)
    plot(t,w(:,1))
    title('Position')
    xlabel('t')
    ylabel('rad')
    
    subplot(2,1,2)
    plot(t,w(:,2))
    title('Velocity')
    xlabel('t')
    ylabel('rad/s')
    
    figure
    plot(w(:,1),w(:,2))
    title('Phase space')
    xlabel('theta')
    ylabel('thetad')
end
end


function rdot = proj(t,r,g,R)
    rdot = [r(2); -g/R*sin(r(1))];
end
