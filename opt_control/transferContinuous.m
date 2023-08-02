%%%%%%%%%%%%%%%%%%%%%
% Continuous function
% Input:
%   Time: input.phase.time - t
%   State: imput.phase.state - x
%       1. Position - r(km): x(i,1:3),3D
%       2. Velocity - v(km/s): x(i,4:6),3D
%       3. Mass - m(kg) : x(i,7),1D
%   Control: input.phase.control - u(kN),3D
%%%%%%%%%%%%%%%%%%%%%
function output=transferContinuous(input)
[~,mu,~]=get_constant();
% Additional constants
Isp=1500;                                   % Specific impulse (s)
g0=9.8e-3;                                  % Gravity constant (km/s2)


% Inputs
t=input.phase.time;
x=input.phase.state;
u=input.phase.control;

% r,v,m, convenient for computing
r=x(:,1:3);
v=x(:,4:6);
m=x(:,7);

% Dynamics
dr=v;
r_norm=norm(r,2);
dv=-mu/r_norm*r+u/m;
u_norm=norm(u,2);
dm=-u_norm/(Isp*g0);
output.phase.dynamics=[dr,dv,dm];

% Path constraint
path=[sum(u.*u,2)];
output.path=path;

% Integrant
output.integrant=norm(u);
end