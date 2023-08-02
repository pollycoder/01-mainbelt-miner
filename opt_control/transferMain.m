clear all;clc

[au,~,day]=get_constant();
[earth_pos,earth_vel]=ephemeris('EARTH');
[mars_pos,mars_vel]=ephemeris('MARS');

X(1)=0.1;
X(2)=1.5;


%%%%%%%%%%%%%%%%%
% Provide data
%%%%%%%%%%%%%%%%%
t0=X(1)*day*365;                        % Initial time (s)
tf=X(2)*day*365;                        % Final time (s)
m0=1.25e3;                              % Initial mass (kg)
rmin=au;                                % Astronic unit (km) - Sun-Earth
rmax=1.6*au;                            % Sun-Mars (km)

% Initial state
index0=round(t0/day);
x0_pos=earth_pos(:,index0);
x0_vel=earth_vel(:,index0);
x0_m=m0;

% Final state
indexf=round(tf/day);
xf_pos=mars_pos(:,indexf);
xf_vel=mars_vel(:,indexf);



%%%%%%%%%%%%%%%
% Boundaries
%%%%%%%%%%%%%%%
% Time bounds
bounds.phase.initialtime.lower=t0;
bounds.phase.initialtime.upper=t0;
bounds.phase.finaltime.lower=tf;
bounds.phase.finaltime.upper=tf;

% State bounds
bounds.phase.initialstate.lower=x0;
bounds.phase.initialstate.upper=x0;


%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%
% Name
setup.name='PlanetTransfer-Problem';                

% Functions
setup.functions.continuous=@transferContinuous;
setup.functions.endpoint=@transferEndpoint;

% Boundaries
setup.bounds=bounds;









