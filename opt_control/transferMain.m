clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer solver - new unit
%   Time: Earth year
%   Distance: AU
%   Mass: Maximum mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------Preparations------------------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[auxdata.au,auxdata.mu,~,auxdata.year]=get_constant(); % (m,m3/s2,s,s)

X(1)=0;
X(2)=20;
Xmid=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for boundaries and normalization 
% Unit transfer:
% mass unit: 
% time unit: Earth year
% length unit: Radius of Earth belt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time
t0=X(1)*auxdata.year;                            % Initial time (s)
tf=X(2)*auxdata.year;                            % Final time (s)
       

% Mass
m0=3e4;                                        % Initial mass (kg)
mf=0.1*m0;                                     % Final minimum residual mass (kg)

% Oribit elements of two orbits
% All the two orbits are circular orbit by default
% Initial orbit elements - Earth
a0=1*auxdata.au;                              % Semi-major axis (m)
e0=0;                                         % Eccentricity (-) - circular
i0=0;                                         % Inclination - plane
Omega0=0;                                     % Longtitude of ascending node (rad) - circular
omega0=0;                                     % Argument of periapsis (rad) - circular
f0=0;                                         % Initial true anomaly (rad) - initial position
[p0,f0,g_0,h0,k0,L0]=coe2mee(a0,e0,i0,Omega0,omega0,f0);

% Interim orbit elements - Mars
% Final orbit elements- Main Belt
aM=1.52*auxdata.au;                           % Semi-major of Earth's orbit (m)
eM=0;                                         % Eccentricity (-) - circular
iM=0;                                         % Inclination (rad) - plane
OmegaM=0;                                     % Longtitude of ascending node (rad) - circular
omegaM=0;                                     % Argument of periapsis (rad) - circular
fM=1.7*pi;                                    % Final true anomaly (rad) - final position
[pM,fM,gM,hM,kM,LM]=coe2mee(aM,eM,iM,OmegaM,omegaM,fM);

% Final orbit elements- Main Belt
af=2.27*auxdata.au;                           % Semi-major of Earth's orbit (m)
ef=0;                                         % Eccentricity (-) - circular
i_f=0;                                        % Inclination (rad) - plane
Omegaf=0;                                     % Longtitude of ascending node (rad) - circular
omegaf=0;                                     % Argument of periapsis (rad) - circular
ff=pi;                                        % Final true anomaly (rad) - final position
[pf,ff,gf,hf,kf,Lf]=coe2mee(af,ef,i_f,Omegaf,omegaf,ff);


%%%%%%%%%%%%%%%%%%
% Unit transfer
%%%%%%%%%%%%%%%%%%
auxdata.m2newUnit=1/auxdata.au;
auxdata.s2newUnit=1/auxdata.year;
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
auxdata.rMarsSOI_newUnit=0.578e9*auxdata.lUnit;
auxdata.muMars_newUnit=4.285e13*auxdata.muUnit;
auxdata.mu_newUnit=auxdata.mu*auxdata.muUnit;                      
auxdata.Tmax_newUnit=[2,2]*auxdata.forceUnit;  
auxdata.g0_newUnit=9.8*auxdata.aUnit;                               
auxdata.Isp_newUnit=3100*auxdata.tUnit;      

t0_newUnit=t0*auxdata.tUnit;
tf_newUnit=tf*auxdata.tUnit;
tMid_newUnit=Xmid;
                                  
% Orbit elements (Spring equinox)
p0_newUnit=a0*auxdata.lUnit;
pM_newUnit=pM*auxdata.lUnit;
pf_newUnit=pf*auxdata.lUnit;
mee0_newUnit=[p0_newUnit,f0,g_0,h0,k0,L0];
meeM_newUnit=[pM_newUnit,fM,gM,hM,kM,LM];
meef_newUnit=[pf_newUnit,ff,gf,hf,kf,Lf];

% Bounds - between Main Belt and Earth
m0_newUnit=m0*auxdata.mUnit;                           
mf_newUnit=mf*auxdata.mUnit;                            
rmin_newUnit=-max(pf,p0)*auxdata.lUnit;                     
rmax_newUnit=-rmin_newUnit;                                 
vmin_newUnit=-30000*auxdata.vUnit;                            
vmax_newUnit=-vmin_newUnit;                                   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------Start solving----------------------------- %                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first phase: from Earth to Mars - all new unit
%
% Escape from Earth and go to Mars
% After calculation, we found that if we escape in the escaping velocity,
% the escaping time will be extremely short, so we just take the Earth
% velocity as our departure velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iphase=1;

% Initial state
x0_m_newUnit=m0_newUnit;  
[x0_pos_newUnit,x0_vel_newUnit]=mee2rv(mee0_newUnit,auxdata.mu_newUnit);
x0_newUnit=[x0_pos_newUnit',x0_vel_newUnit',x0_m_newUnit];    

% Interim state
mM_newUnit=0.85;               % Still have 70% of mass left          
xM_m_newUnit=mM_newUnit;

[Mars_pos_newUnit,Mars_vel_newUnit]=mee2rv(meeM_newUnit,auxdata.mu_newUnit);    % Mars position and velocity

% GA-Velocity
% Direction
MarsVelNorm=sqrt(Mars_vel_newUnit(1)^2+Mars_vel_newUnit(2)^2+Mars_vel_newUnit(3)^2);
unitVec=Mars_vel_newUnit./MarsVelNorm;

% Position
rRelative_newUnit=-auxdata.rMarsSOI_newUnit*unitVec;
Mars_GA_pos_newUnit=rRelative_newUnit+Mars_pos_newUnit;

% Amount and direction of velocity
rotate=[cos(pi/3),-sin(pi/3),0;sin(pi/3),cos(pi/3),0;0,0,1];
vRelativeNorm_newUnit=sqrt(4/sqrt(3)*auxdata.muMars_newUnit/auxdata.rMarsSOI_newUnit);
unitVec=rotate*unitVec;

% Final velocity
Mars_GA_vel_newUnit=Mars_vel_newUnit+vRelativeNorm_newUnit*unitVec;
auxdata.MarsState_newUnit=[Mars_pos_newUnit',Mars_vel_newUnit'];
  
xM_newUnit=[Mars_GA_pos_newUnit',Mars_GA_vel_newUnit',xM_m_newUnit];

% Guess
NGuess(1)=1;                            % Number of loops
N=NGuess(1);
Node=100;                               % Nodes of one loop
guessNode=Node*(N+1);                   % Total guess nodes

meeM_newUnit(6)=mod(meeM_newUnit(6),2*pi)+N*2*pi;

for i=1:guessNode
    meetemp_newUnit=mee0_newUnit-(mee0_newUnit-meeM_newUnit)/(guessNode-1)*(i-1);
    meetemp_newUnit(1)=abs(meetemp_newUnit(1));
    meetemp_newUnit(2)=abs(meetemp_newUnit(2));
    [rtemp_newUnit,vtemp_newUnit]=mee2rv(meetemp_newUnit,auxdata.mu_newUnit);
    mtemp_newUnit=m0_newUnit-(m0_newUnit-mM_newUnit)/(guessNode-1)*(i-1);
    guess.phase(iphase).state(i,:)=[rtemp_newUnit',vtemp_newUnit',mtemp_newUnit];
    guess.phase(iphase).control(i,:)=[0.1,0.1,0];
end
guess.phase(iphase).time=(t0_newUnit:(tMid_newUnit-t0_newUnit)/(guessNode-1):tMid_newUnit)';

% Time bounds
bounds.phase(iphase).initialtime.lower=t0_newUnit;
bounds.phase(iphase).initialtime.upper=t0_newUnit;
bounds.phase(iphase).finaltime.lower=t0_newUnit;
bounds.phase(iphase).finaltime.upper=tMid_newUnit;

% State bounds
bounds.phase(iphase).initialstate.lower=x0_newUnit;
bounds.phase(iphase).initialstate.upper=x0_newUnit;
bounds.phase(iphase).finalstate.lower=[xM_newUnit(1:6),xM_newUnit(7)];
bounds.phase(iphase).finalstate.upper=[xM_newUnit(1:6),x0_newUnit(7)];
bounds.phase(iphase).state.lower=[rmin_newUnit*ones(1,3),vmin_newUnit*ones(1,2),0,mf_newUnit];
bounds.phase(iphase).state.upper=[rmax_newUnit*ones(1,3),vmax_newUnit*ones(1,2),0,m0_newUnit];

% Control bounds
bounds.phase(iphase).control.lower=-[ones(1,2),0];
bounds.phase(iphase).control.upper=+[ones(1,2),0];

% Path constraint
bounds.phase(iphase).path.lower=0;
bounds.phase(iphase).path.upper=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second phase: from Mars to Main Belt - all new unit
% Now we go to the Main Belt. 
% We picked one fixed position to stop.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iphase=2;

% Initial time
t02_newUnit=tMid_newUnit;
tf2_newUnit=tf_newUnit;

% Initial state - only mass
x02_m_newUnit=mM_newUnit;

% Final state          
xf_m_newUnit=mf_newUnit;
[xf_pos_newUnit,xf_vel_newUnit]=mee2rv(meef_newUnit,auxdata.mu);
xf_newUnit=[xf_pos_newUnit',xf_vel_newUnit',xf_m_newUnit];

% Guess
NGuess(2)=10;                            % Number of loops
N=NGuess(2);
Node=100;                               % Nodes of one loop
guessNode=Node*(N+1);                   % Total guess nodes

meeM_newUnit(6)=mod(meeM_newUnit(6),2*pi);
meef_newUnit(6)=mod(meef_newUnit(6),2*pi)+N*2*pi;

for i=1:guessNode
    meetemp_newUnit=meeM_newUnit-(meeM_newUnit-meef_newUnit)/(guessNode-1)*(i-1);
    meetemp_newUnit(1)=abs(meetemp_newUnit(1));
    meetemp_newUnit(2)=abs(meetemp_newUnit(2));
    [rtemp_newUnit,vtemp_newUnit]=mee2rv(meetemp_newUnit,auxdata.mu_newUnit);
    mtemp_newUnit=mM_newUnit-(mM_newUnit-mf_newUnit)/(guessNode-1)*(i-1);
    guess.phase(iphase).state(i,:)=[rtemp_newUnit',vtemp_newUnit',mtemp_newUnit];
    guess.phase(iphase).control(i,:)=[0.1,0.1,0];
end
guess.phase(iphase).time=(t02_newUnit:(tf_newUnit-t02_newUnit)/(guessNode-1):tf_newUnit)';

% Time bounds
bounds.phase(iphase).initialtime.lower=t0_newUnit;
bounds.phase(iphase).initialtime.upper=tf_newUnit;
bounds.phase(iphase).finaltime.lower=t0_newUnit;
bounds.phase(iphase).finaltime.upper=tf_newUnit;

% State bounds
bounds.phase(iphase).initialstate.lower=[xM_newUnit(1:3),vmin_newUnit*ones(1,2),0,xM_newUnit(7)];
bounds.phase(iphase).initialstate.upper=[xM_newUnit(1:3),vmax_newUnit*ones(1,2),0,xM_newUnit(7)];
bounds.phase(iphase).finalstate.lower=[xf_newUnit(1:6),xf_newUnit(7)];
bounds.phase(iphase).finalstate.upper=[xf_newUnit(1:6),xM_newUnit(7)];
bounds.phase(iphase).state.lower=[rmin_newUnit*ones(1,3),vmin_newUnit*ones(1,2),0,mf_newUnit];
bounds.phase(iphase).state.upper=[rmax_newUnit*ones(1,3),vmax_newUnit*ones(1,2),0,mM_newUnit];

% Control bounds
bounds.phase(iphase).control.lower=-[ones(1,2),0];
bounds.phase(iphase).control.upper=+[ones(1,2),0];

% Path constraint
bounds.phase(iphase).path.lower=0;
bounds.phase(iphase).path.upper=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bounds for events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bounds.eventgroup.lower=zeros(1,8);
bounds.eventgroup.upper=zeros(1,8);

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
setup.mesh.colpointsmax = 1000;
setup.mesh.colpointsmin = 100;
for i=1:iphase
setup.mesh.phase(i).colpoints = 10*ones(1,10);
setup.mesh.phase(i).fraction =  0.1*ones(1,10);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Solution - start gpops
%%%%%%%%%%%%%%%%%%%%%%%%
output=gpops2(setup);

result=output.result.objective;
fprintf("J=%f",result);


% Plot
u1=output.result.solution.phase(1).control;
u2=output.result.solution.phase(2).control;
u=[u1;u2];
ux=u2(:,1);
uy=u2(:,2);
uz=u2(:,3);
unorm=sqrt(u2(:,1).^2+u2(:,2).^2+u2(:,3).^2);
r1=output.result.solution.phase(1).state(:,1:3);
r2=output.result.solution.phase(2).state(:,1:3);
r=[r1;r2];
x=r2(:,1);
y=r2(:,2);
z=r2(:,3);
x1=r2(:,1);
y1=r2(:,2);
z1=r2(:,3);
t1=output.result.solution.phase(1).time;
t2=output.result.solution.phase(2).time;
t=[t1,t2];

figure
quiver3(x,y,z,ux,uy,uz,0.3,'r');hold on;
plot3(x,y,z,'LineWidth',1.5,'Color','black');hold on;
plot3(x(1),y(1),z(1),'r*','LineWidth',2);hold on;
text(x(1),y(1),z(1),'Departure');hold on;
plot3(x(end),y(end),z(end),'g*','LineWidth',2);hold on;
text(x(end),y(end),z(end),'Arrival');hold on;
plot3(x1(1),y1(1),z1(1),'r*','LineWidth',2);hold on;
text(x1(1),y1(1),z1(1),'GA-point');hold on;

axis equal

figure
plot(t2,unorm,'LineWidth',1.5,'Color','b');


