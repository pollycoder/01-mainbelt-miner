%%%%%%%%%%%%%%%%%%%%%
% Continuous function - all new unit
% Input:
%   Time: input.phase.time - t
%   State: imput.phase.state - x
%       1. Position - r: x(i,1:3),3D
%       2. Velocity - v: x(i,4:6),3D
%       3. Mass - m : x(i,7),1D
%   Control: input.phase.control - u,3D
%%%%%%%%%%%%%%%%%%%%%
function output=transferContinuous(input)
% Additional constants
Isp_newUnit=input.auxdata.Isp_newUnit;                         
g0_newUnit=input.auxdata.g0_newUnit;       
Tmax_newUnit=input.auxdata.Tmax_newUnit;                            


% Inputs
x_newUnit=input.phase.state;
u_newUnit=input.phase.control;                      

% r,v,m, convenient for computing
r_newUnit=x_newUnit(:,1:3);
v_newUnit=x_newUnit(:,4:6);
m_newUnit=x_newUnit(:,7);

% Dynamics
dr_newUnit=v_newUnit;          

num_points=size(r_newUnit,1);
acc1_newUnit=zeros(num_points,3);


for i=1:num_points
    rt_newUnit=r_newUnit(i,:);
    vt_newUnit=v_newUnit(i,:);

    % Transfer to international unit
    rt=rt_newUnit./input.auxdata.lUnit; % (m)
    vt=vt_newUnit./input.auxdata.vUnit; % (m/s)
    rt_norm=sqrt(rt(1).^2+rt(2).^2+rt(3).^2); % (m)

    acc1=-input.auxdata.mu./rt_norm.^3.*rt;
    
    % Transfer back to au and day
    acc1_newUnit(i,:)=acc1.*input.auxdata.aUnit;
end

acc2_newUnit=diag(m_newUnit.^-1)*u_newUnit.*Tmax_newUnit;                

dv_newUnit=acc1_newUnit+acc2_newUnit;

unorm_newUnit=sqrt(u_newUnit(:, 1).^2 + u_newUnit(:, 2).^2 + u_newUnit(:, 3).^2);
dm_newUnit=-Tmax_newUnit.*unorm_newUnit./(Isp_newUnit.*g0_newUnit);
output.dynamics=[dr_newUnit,dv_newUnit,dm_newUnit];

% Path constraint
output.path=unorm_newUnit;
end