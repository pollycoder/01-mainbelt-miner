%%%%%%%%%%%%%%%%%%%
% Two-body problem
%%%%%%%%%%%%%%%%%%%
range=[0,1;0,5];


[X,result]=pso(@(x)obj_func(x),range,10000,1000);


% X(0):start time
% X(1):end time
function result=obj_func(X)
earth_T=24*3600*365;                    % Earth year (s)
mars_T=24*3600*686.971;                 % Mar year (s)

delta_X=X(2)-X(1);                      % Normalized total_time
start_time=X(1)*earth_T;                % Waiting time
total_time=delta_X*earth_T;

au=149597582503.1;                      % 1AU(m)
earth_R=au;                             % Radius of Earth-Sun orbit(m)
mars_R=1.52*au;                         % Radius of Mars-Sun orbit(m)

earth_angv=2*pi/earth_T;                % Angular velocity of Earth
mars_angv=2*pi/mars_T;                  % Angular velocity of mars

% Initial position
earth_theta=pi;
mars_theta=2/3*pi;
earth_theta=earth_theta+start_time*earth_angv;
earth_position=earth_R*[cos(earth_theta),sin(earth_theta),0];

% Final position
mars_theta=mars_theta+total_time*mars_angv;
mars_position=earth_R*[cos(mars_theta),sin(mars_theta),0];

% Lambert
[v1,v2]=LAMBERTBATTIN(earth_position,mars_position,'pro',total_time);

v0=earth_angv*earth_R*[sin(earth_theta),cos(earth_theta),0];
vt=mars_angv*mars_R*[sin(mars_theta),cos(mars_theta),0];
dv1=delta_v(v1,v0);
dv2=delta_v(v2,vt);
result=dv1+dv2;
end
