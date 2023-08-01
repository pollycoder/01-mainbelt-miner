%%%%%%%%%%%%%%%%%%%%%%%%
% Object function
%%%%%%%%%%%%%%%%%%%%%%%%

% X(0):start time
% X(1):end time
function result=obj_func(X)
[~,mu_sun,day]=get_constant();

% Set the time
delta_X=X(2)-X(1);                      % Normalized transfer_time
if delta_X<=0 || delta_X>5              % Penalty function
    result=1E20;
    return
end
wait_time=X(1)*day*365;                 % Waiting time (s)
transfer_time=delta_X*day*365;          % Transfer time (s)
total_time=X(2)*day*365;                % Total time (s)

% Read the ephemeris
[earth_pos,earth_vel]=ephemeris('EARTH');
[mars_pos,mars_vel]=ephemeris('MARS');

% Initial position
index=round(wait_time/day)+1;
r0=earth_pos(:,index);
v0=earth_vel(:,index);

% Final position
index=round(total_time/day)+1;
rf=mars_pos(:,index);
vf=mars_vel(:,index);

% Lambert
[v1,v2,~,~,~,~]=LambSol(r0,rf,transfer_time,mu_sun);

dv1=delta_v(v1,v0);
dv2=delta_v(v2,vf);
result=dv1+dv2;
end
