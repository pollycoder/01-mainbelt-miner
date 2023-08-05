clear all;clc

[auxdata.au,auxdata.mu,auxdata.day]=get_constant();
%[init_pos,init_vel]=ephemeris('EARTH')*1000;
%[final_pos,final_vel]=ephemeris('MARS')*1000;

phi0=0;dphi=8.2;
Rbelt=2.27*auxdata.au;                              % (m)
REarth=1*auxdata.au;                                % (m)
Vbelt=sqrt(auxdata.mu/Rbelt);                       % (m/s)
VEarth=sqrt(auxdata.mu/REarth);                     % (m/s)
x0_pos=Rbelt*[cos(phi0);sin(phi0);0];               % (m)
x0_vel=Vbelt*[-sin(phi0);cos(phi0);0];              % (m/s)
xf_pos=REarth*[cos(phi0+dphi);sin(phi0+dphi);0];    % (m)
xf_vel=VEarth*[-sin(phi0+dphi);cos(phi0+dphi);0];   % (m/s)    

X(1)=0;
X(2)=5;


%%%%%%%%%%%%%%%%%
% Provide data
%%%%%%%%%%%%%%%%%
t0=X(1)*365*auxdata.day;                % Initial time (s)
tf=X(2)*365*auxdata.day;                % Final time (s)
m0=10000000.0000;                       % Initial mass (kg)
mf=700.0000;                            % Final minimum mass (kg)
rmin=-2.27*auxdata.au;                  % Minimum R (m)
rmax=-rmin;                             % Maximum R (m)
vmin=-40000;                            % Maximum V (m/s)
vmax=-vmin;                             % Minimum V (m/s)

% Initial state
%x0_pos=init_pos(:,round(t0)+1)./auxdata.au;
%x0_vel=2*init_vel(:,round(t0)+1)./auxdata.au.*auxdata.day;
x0_m=m0;
x0=[x0_pos',x0_vel',x0_m];

% Final state
indexf=round(tf)+1;
%xf_pos=final_pos(:,indexf)./auxdata.au;
%xf_vel=final_vel(:,indexf)./auxdata.au.*auxdata.day;
xf_m=mf;
xf=[xf_pos',xf_vel',xf_m];


%%%%%%%%%%%%%%%%%%%%%%%
% Auxdata
%%%%%%%%%%%%%%%%%%%%%%%
auxdata.Tmax=300;                                   % Maximum thrust (kg*m/s2)
auxdata.g0=9.8;                                     % Gravity coefficient (m/s2)
auxdata.Isp=1000;                                   % Specific impulse (s)



%%%%%%%%%%%%%%%
% Boundaries
%%%%%%%%%%%%%%%
% Time bounds
bounds.phase.initialtime.lower=t0;
bounds.phase.initialtime.upper=t0;
bounds.phase.finaltime.lower=t0;
bounds.phase.finaltime.upper=tf;

% State bounds
bounds.phase.initialstate.lower=[x0(1:3),vmin*ones(1,3),m0];
bounds.phase.initialstate.upper=[x0(1:3),vmax*ones(1,3),m0];
bounds.phase.finalstate.lower=[xf(1:3),vmin*ones(1,3),mf];
bounds.phase.finalstate.upper=[xf(1:3),vmax*ones(1,3),m0];
bounds.phase.state.lower=[rmin*ones(1,3),vmin*ones(1,3),mf];
bounds.phase.state.upper=[rmax*ones(1,3),vmax*ones(1,3),m0];

% Control bounds
bounds.phase.control.lower=-ones(1,3);
bounds.phase.control.upper=+ones(1,3);

% Path constraint
bounds.phase.path.lower=0;
bounds.phase.path.upper=1;

% Integral constraint
bounds.phase.integral.lower=t0;
bounds.phase.integral.upper=tf;


%%%%%%%%%%%%%%%
% Guess
%%%%%%%%%%%%%%%
% Time guess
guess.phase.time=[t0;tf];

% State guess
guess.phase.state(:,1)=[x0_pos(1);x0_pos(1)];
guess.phase.state(:,2)=[x0_pos(2);x0_pos(2)];
guess.phase.state(:,3)=[x0_pos(3);x0_pos(3)];
guess.phase.state(:,4)=[x0_vel(1);x0_vel(1)];
guess.phase.state(:,5)=[x0_vel(2);x0_vel(2)];
guess.phase.state(:,6)=[x0_vel(3);x0_vel(3)];
guess.phase.state(:,7)=[x0_m;x0_m];

% Control guess
guess.phase.control(:,1)=[auxdata.Tmax;auxdata.Tmax];
guess.phase.control(:,2:4)=-[x0_vel';x0_vel']/norm(x0_vel,2);

% Integral guess
guess.phase.integral=0;

%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%
% Name
setup.name='LowThrustTransfer-Problem';                

% Functions
setup.functions.continuous=@transferContinuous;
setup.functions.endpoint=@transferEndpoint;

% Auxdatas
setup.auxdata=auxdata;

% Boundaries
setup.bounds=bounds;

% Guess
setup.guess=guess;

% Solver
setup.nlp.solver='snopt';

% Derivative
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';

% Mesh
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-5;
setup.mesh.maxiteration = 5;
setup.mesh.colpointsmax = 4;
setup.mesh.colpointsmin = 10;
setup.mesh.phase.colpoints = 4*ones(1,10);
setup.mesh.phase.fraction =  0.1*ones(1,10);

%%%%%%%%%%%%%%%%%%%%%%%%
% Solution - start gpops
%%%%%%%%%%%%%%%%%%%%%%%%
output=gpops2(setup);
result=output.result.objective;
fprintf("J=%f",result);

x=output.result.solution.phase.state(:,1);
y=output.result.solution.phase.state(:,2);
u=output.result.solution.phase.control(:,1);
t=output.result.solution.phase.time;
figure
plot(x,y);
figure
plot(t,u);