%%%%%%%%%%%%%%%%%%%%%
% Continuous function
% Input:
%   Time: input.phase.time - t
%   State: imput.phase.state - x
%       1. Position - r(m): x(i,1:3),3D
%       2. Velocity - v(m/s): x(i,4:6),3D
%       3. Mass - m(kg) : x(i,7),1D
%   Control: input.phase.control - u(N),3D
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
r_norm=sqrt(r(:,1).^2+r(:,2).^3+r(:,3)^2);                  % ||r||

acc1=-input.auxdata.mu./(r_norm.^3).*r;                     % Gravity
acc2=diag(1./m)*u.*Tmax;                                    % Thrust
dv=acc1+acc2;
unorm = sqrt(u(:, 1).^2 + u(:, 2).^2 + u(:, 3).^2);
dm=-Tmax.*unorm./(Isp.*g0);
output.dynamics=[dr,dv,dm];

% Path constraint
output.path=unorm;
output.integrand=Tmax./(Isp.*g0).*unorm;
end