clear all;clc

[auxdata.au,auxdata.mu,auxdata.day]=get_constant(); % (m,m3/s2,s)

X(1)=0;
X(2)=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for boundaries and normalization 
% Unit transfer:
% mass unit: m0
% time unit: Earth year
% length unit: Radius of Earth belt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time
t0=X(1)*365*auxdata.day;                % Initial time (s)
tf=X(2)*365*auxdata.day;                % Final time (s)
       

% Mass
m0=1e7;                                 % Initial mass (kg)
mf=0.3*m0;                              % Final minimum residual mass (kg)

% Oribit elements of two orbits
% All the two orbits are circular orbit by default
% Main belt
aBelt=2.27*auxdata.au;                  % Semi-major axis of main belt (m)
eBelt=0;                                % Eccentricity of main belt (-) - circular
iBelt=0;                                % Inclination of main belt (rad) - plane
OmegaBelt=0;                            % Longtitude of ascending node (rad) - circular
omegaBelt=0;                            % Argument of periapsis (rad) - circular
fBelt_0=0;                              % Initial true anomaly (rad) - initial position

% Earth
aEarth=1*auxdata.au;                    % Semi-major of Earth's orbit (m)
eEarth=0;                               % Eccentricity s of Earth's orbit (-) - circular
iEarth=0;                               % Inclination of Earth's orbit (rad) - plane
OmegaEarth=0;                           % Longtitude of Earth's orbit (rad) - circular
omegaEarth=0;                           % Argument of periapsis (rad) - circular
fEarth_f=8.2-2*pi;                      % Final true anomaly (rad) - final position
TOrbit=2*pi*sqrt(aEarth^3/auxdata.mu);  % Earth year (s)


%%%%%%%%%%%%%%%%%%
% Unit transfer
%%%%%%%%%%%%%%%%%%
auxdata.m2newUnit=1/aEarth;
auxdata.s2newUnit=1/(TOrbit);
auxdata.tUnit=auxdata.s2newUnit;
auxdata.lUnit=auxdata.m2newUnit;  
auxdata.mUnit=1/m0;
auxdata.vUnit=auxdata.m2newUnit/auxdata.s2newUnit;
auxdata.aUnit=auxdata.m2newUnit/auxdata.s2newUnit^2;
auxdata.muUnit=auxdata.m2newUnit^3/auxdata.s2newUnit^2;
auxdata.forceUnit=auxdata.aUnit*auxdata.mUnit;


%%%%%%%%%%%%%%%%%
% Provide data - all new unit
%%%%%%%%%%%%%%%%%
auxdata.mu_newUnit=auxdata.mu*auxdata.muUnit;                      
auxdata.Tmax_newUnit=2000*auxdata.forceUnit;                        
auxdata.g0_newUnit=9.8*auxdata.aUnit;                               
auxdata.Isp_newUnit=3100*auxdata.tUnit;      

t0_newUnit=t0*auxdata.tUnit;
tf_newUnit=tf*auxdata.tUnit;
                                  
% Orbit elements
aBelt_newUnit=aBelt*auxdata.lUnit;
aEarth_newUnit=aEarth*auxdata.lUnit;
coe0_newUnit=[aBelt_newUnit,eBelt,iBelt,OmegaBelt,omegaBelt,fBelt_0];
coef_newUnit=[aEarth_newUnit,eEarth,iEarth,omegaEarth,OmegaEarth,fEarth_f];

% Bounds
m0_newUnit=m0*auxdata.mUnit;                           
mf_newUnit=mf*auxdata.mUnit;                            
rmin_newUnit=-2.27*auxdata.au*auxdata.lUnit;                     
rmax_newUnit=-rmin_newUnit;                                 
vmin_newUnit=-30000*auxdata.vUnit;                            
vmax_newUnit=-vmin_newUnit;                                   

% Initial state
x0_m_newUnit=m0_newUnit;  
[x0_pos_newUnit,x0_vel_newUnit]=coe2rv(coe0_newUnit,auxdata.mu_newUnit);
x0_newUnit=[x0_pos_newUnit',x0_vel_newUnit',x0_m_newUnit];    

% Final state
xf_m_newUnit=mf_newUnit; 
[xf_pos_newUnit,xf_vel_newUnit]=coe2rv(coef_newUnit,auxdata.mu_newUnit);
xf_newUnit=[xf_pos_newUnit',xf_vel_newUnit',xf_m_newUnit];    


%%%%%%%%%%%%%%%%%%%%
% Guess
%%%%%%%%%%%%%%%%%%%%
NGuess=2;                               % Number of loops

N=NGuess;
Node=1000;                              % Nodes of one loop
guessNode=Node*(N+1);                   % Total guess nodes

coef_newUnit(6)=mod(coef_newUnit(6),2*pi)+N*2*pi;

for i=1:guessNode
    coetemp_newUnit=coe0_newUnit-(coe0_newUnit-coef_newUnit)/(guessNode-1)*(i-1);
    coetemp_newUnit(1)=abs(coetemp_newUnit(1));
    coetemp_newUnit(2)=abs(coetemp_newUnit(2));
    [rtemp_newUnit,vtemp_newUnit]=coe2rv(coetemp_newUnit,auxdata.mu_newUnit);
    mtemp_newUnit=m0_newUnit-(m0_newUnit-mf_newUnit)/(guessNode-1)*(i-1);
    guess.phase.state(i,:)=[rtemp_newUnit',vtemp_newUnit',mtemp_newUnit];
    guess.phase.control(i,:)=[0.1,0.1,0];
end
guess.phase.time=(t0_newUnit:(tf_newUnit-t0_newUnit)/(guessNode-1):tf_newUnit)';


%%%%%%%%%%%%%%%
% Boundaries - all new unit
%%%%%%%%%%%%%%%
% Time bounds
bounds.phase.initialtime.lower=t0_newUnit;
bounds.phase.initialtime.upper=t0_newUnit;
bounds.phase.finaltime.lower=t0_newUnit;
bounds.phase.finaltime.upper=tf_newUnit;

% State bounds
bounds.phase.initialstate.lower=x0_newUnit;
bounds.phase.initialstate.upper=x0_newUnit;
bounds.phase.finalstate.lower=[xf_newUnit(1:6),xf_newUnit(7)];
bounds.phase.finalstate.upper=[xf_newUnit(1:6),x0_newUnit(7)];
bounds.phase.state.lower=[rmin_newUnit*ones(1,3),vmin_newUnit*ones(1,2),0,mf_newUnit];
bounds.phase.state.upper=[rmax_newUnit*ones(1,3),vmax_newUnit*ones(1,2),0,m0_newUnit];

% Control bounds
bounds.phase.control.lower=-[ones(1,2),0];
bounds.phase.control.upper=+[ones(1,2),0];

% Path constraint
bounds.phase.path.lower=0;
bounds.phase.path.upper=1;


%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%
% Name
setup.name='LowThrustTransfer-Problem'; 
setup.method='RPMintegration';
setup.scale.method='automatic-bounds';
setup.tolerance=1e-9;

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
setup.nlp.snoptoptions.maxiterations=1e6;

% Derivative
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';

% Mesh
setup.mesh.method = 'hp-LiuRao-Legendre';
setup.mesh.tolerance = 1e-9;
setup.mesh.maxiteration = 10;
setup.mesh.colpointsmax = 10000;
setup.mesh.colpointsmin = 1000;
setup.mesh.phase.colpoints = 10*ones(1,10);
setup.mesh.phase.fraction =  0.1*ones(1,10);

%%%%%%%%%%%%%%%%%%%%%%%%
% Solution - start gpops
%%%%%%%%%%%%%%%%%%%%%%%%
output=gpops2(setup);
result=output.result.objective;
fprintf("J=%f",result);

% Plot
u=output.result.solution.phase.control;
ux=u(:,1);
uy=u(:,2);
uz=u(:,3);
unorm=sqrt(u(:,1).^2+u(:,2).^2+u(:,3).^2);
r=output.result.solution.phase.state(:,1:3);
x=r(:,1);
y=r(:,2);
z=r(:,3);
t=output.result.solution.phase.time;

subplot(1,2,1);
quiver3(x,y,z,ux,uy,uz,2.5,'r');hold on;
plot3(x,y,z,'LineWidth',1.5,'Color','black');hold on;
plot3(x(1),y(1),z(1),'r*','LineWidth',2);hold on;
text(x(1),y(1),z(1),'Departure');hold on;
plot3(x(end),y(end),z(end),'g*','LineWidth',2);hold on;
text(x(end),y(end),z(end),'Arrival');
axis equal

subplot(1,2,2);
plot(t,unorm,'LineWidth',1.5,'Color','b');