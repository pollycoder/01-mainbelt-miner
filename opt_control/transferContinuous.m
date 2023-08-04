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
% Additional constants
Isp=input.auxdata.Isp;                         
g0=input.auxdata.g0;       
Tmax=input.auxdata.Tmax;                            


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
r_norm=sqrt(sum(r.*r,2));
dv=-input.auxdata.mu./(r_norm.^3).*r+u(:,1).*u(:,2:4)./m;
dm=-u(:,1)./(Isp*g0);
output.dynamics=[dr,dv,dm];

% Path constraint
path=[sum(u(:,2:4).*u(:,2:4),2)];
output.path=path;
output.integrand=sqrt(sum(u(:,1).*u(:,1),2));
end