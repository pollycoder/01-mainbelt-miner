clear all
clc
%% Transfer solver - new unit
%   Time: Earth year
%   Distance: AU
%   Mass: Maximum mass




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------Preparations------------------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic elements
[auxdata.au,auxdata.mu,~,auxdata.year]=get_constant(); % (m,m3/s2,s,s)t
% Time points (y)
t0=0;
tf=200;
tMid1=4;
tMid2=5;

% Mass requirements
m0=3e4;                                                 % Initial mass (kg)

% Parameter of the spacecraft
Isp=4000;                                               % Specific impulse (s)
Tmax=2;                                                 % Maximum thrust (N)
g0=9.8;                                                 % Gravity acceleration (m/s2)

% Parameter of the planets
% Mars
muMars=6.67e-11*6.4219e23;                              % Gravity coefficient of Mars (m3/s2)
rMarsSOI=0.578e9;                                       % SOI radius of Mars (m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Unit transfer and provide auxillary data
% mass unit: proportion of initial total mass (m0)
% time unit: Earth year
% length unit: AU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxdata.m2newUnit=1/auxdata.au;
auxdata.s2newUnit=1/auxdata.year;

auxdata.tUnit=auxdata.s2newUnit;
auxdata.lUnit=auxdata.m2newUnit;  
auxdata.mUnit=1/m0;
auxdata.vUnit=auxdata.m2newUnit/auxdata.s2newUnit;
auxdata.aUnit=auxdata.m2newUnit/auxdata.s2newUnit^2;
auxdata.muUnit=auxdata.m2newUnit^3/auxdata.s2newUnit^2;
auxdata.forceUnit=auxdata.aUnit*auxdata.mUnit;

auxdata.Isp_newUnit=Isp*auxdata.tUnit;                  % Specific impulse
auxdata.g0_newUnit=g0*auxdata.aUnit;                    % Gravity acceleration
auxdata.Tmax_newUnit=Tmax*auxdata.forceUnit;            % Maximum thrust
auxdata.mu_newUnit=auxdata.mu*auxdata.muUnit;           % Gravity coefficient of Sun
auxdata.muMars_newUnit=muMars*auxdata.muUnit;   % Gravity coeffient of Mars
auxdata.rMarsSOI_newUnit=rMarsSOI*auxdata.lUnit;% SOI radius of Mars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descriptions of the 3 orbits - Earth, Mars, Main Belt
% Use spring squinox orbit elements 
% All in new unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial orbit elements - Earth
a0_newUnit=1;                                 % Semi-major axis (AU)
e0=0.0167;                                         % Eccentricity (-) - circular
i0=0;                                         % Inclination - plane
Omega0=0;                                     % Longtitude of ascending node (rad) - circular
omega0=0;                                     % Argument of periapsis (rad) - circular
f0=0;                                         % Initial true anomaly (rad) - initial position
[p0_newUnit,f0,g_0,h0,k0,L0]=coe2mee(a0_newUnit,e0,i0,Omega0,omega0,f0);
auxdata.mee0_newUnit=[p0_newUnit,f0,g_0,h0,k0,L0];

% Interim orbit elements - Mars
aM_newUnit=1.52;                              % Semi-major of Earth's orbit (AU)
eM=0.0934;                                    % Eccentricity (-) - circular
iM=0.0323;                                    % Inclination (rad) - plane
OmegaM=0;                                     % Longtitude of ascending node (rad) - circular
omegaM=0;                                     % Argument of periapsis (rad) - circular
fM=1.5*pi;                                    % Final true anomaly (rad) - final position
[pM_newUnit,fM,gM,hM,kM,LM]=coe2mee(aM_newUnit,eM,iM,OmegaM,omegaM,fM);
auxdata.meeM_newUnit=[pM_newUnit,fM,gM,hM,kM,LM];

% Final orbit elements- Main Belt
af_newUnit=2.27;                              % Semi-major of Earth's orbit (AU)
ef=0;                                         % Eccentricity (-) - circular
i_f=0;                                        % Inclination (rad) - plane
Omegaf=0;                                     % Longtitude of ascending node (rad) - circular
omegaf=0;                                     % Argument of periapsis (rad) - circular
ff=1.8*pi;                                    % Final true anomaly (rad) - final position
[pf_newUnit,ff,gf,hf,kf,Lf]=coe2mee(af_newUnit,ef,i_f,Omegaf,omegaf,ff);
auxdata.meef_newUnit=[pf_newUnit,ff,gf,hf,kf,Lf];


% Bounds
mMax_newUnit=1;                                 % m0 (unit)
mMid_newUnit=0.8;                               % Mass on Mars (unit),still have >=80% left
mMin_newUnit=0.3;                               % Final mass (unit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------Start solving----------------------------- %                                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iphase=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first phase: from Earth to Mars - all new unit
%
% Escape from Earth and go to Mars
% After calculation, we found that if we escape in the escaping velocity,
% the escaping time will be extremely short, so we just take the Earth
% velocity as our departure velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guess
NGuess=2;                            % Number of loops
N=NGuess(1);
Node=100;                             % Nodes of one loop
guessNode=Node*(N+1);                   % Total guess nodes

auxdata.meeM_newUnit(6)=mod(auxdata.meeM_newUnit(6),2*pi)+N*2*pi;

% Bounds
Lmax=auxdata.meeM_newUnit(6);                   % L (rad)                                        
Lmin=auxdata.mee0_newUnit(6);
fmax=fM;
fmin=0;
pmax_newUnit=auxdata.meeM_newUnit(1);               % p (AU)
pmin_newUnit=auxdata.mee0_newUnit(1)*(1-0.3^2);
amin=a0_newUnit;
amax=aM_newUnit;
hmin=auxdata.mee0_newUnit(4);
hmax=auxdata.meeM_newUnit(4);


% Start and end
Lstart=auxdata.mee0_newUnit(6);
Lend=auxdata.meeM_newUnit(6);
pstart_newUnit=auxdata.mee0_newUnit(1);
pend_newUnit=auxdata.meeM_newUnit(1);
hstart=auxdata.mee0_newUnit(4);
hend=auxdata.meeM_newUnit(4);
mstart_newUnit=mMax_newUnit;
mend_newUnit=mMid_newUnit;

for i=1:guessNode
    ptemp_newUnit=pstart_newUnit+(pend_newUnit-pstart_newUnit)*(i-1)/(guessNode-1);
    Ltemp=Lstart+(Lend-Lstart)*(i-1)/(guessNode-1);
    htemp=hstart+(hend-hstart)*(i-1)/(guessNode-1);
    mtemp_newUnit=mstart_newUnit+(mend_newUnit-mstart_newUnit)*(i-1)/(guessNode-1);
    guess.phase(iphase).state(i,:)=[ptemp_newUnit,0,0,0,0,Ltemp,mtemp_newUnit];
    guess.phase(iphase).control(i,:)=[0.1,0.1,1e-3];
end
guess.phase(iphase).time=(t0:(tMid1-t0)/(guessNode-1):tMid1)';

% Time bounds
bounds.phase(iphase).initialtime.lower=t0;
bounds.phase(iphase).initialtime.upper=t0;
bounds.phase(iphase).finaltime.lower=t0;
bounds.phase(iphase).finaltime.upper=tMid1;

% State bounds
bounds.phase(iphase).initialstate.lower=[pstart_newUnit,f0,0,hstart,0,Lstart,mMax_newUnit];
bounds.phase(iphase).initialstate.upper=[pstart_newUnit,f0,0,hstart,0,Lstart,mMax_newUnit];
bounds.phase(iphase).finalstate.lower=[pend_newUnit,fM,0,hend,0,Lmin,mMid_newUnit];
bounds.phase(iphase).finalstate.upper=[pend_newUnit,fM,0,hend,0,Lmax,mMax_newUnit];
bounds.phase(iphase).state.lower=[pmin_newUnit,fmin,0,hmin,0,Lmin,mMid_newUnit];
bounds.phase(iphase).state.upper=[pmax_newUnit,fmax,0,hmax,0,Lmax,mMax_newUnit];

% Control bounds
bounds.phase(iphase).control.lower=-ones(1,3);
bounds.phase(iphase).control.upper=+ones(1,3);

% Path constraint
bounds.phase(iphase).path.lower=0;
bounds.phase(iphase).path.upper=1;



%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%
% Name
setup.name='LowThrustTransfer-Problem'; 
setup.method='RPMintegration';
setup.scale.method='automatic-bounds';
setup.tolerance=1e-9;

% Functions
setup.functions.continuous=@transferGAContinuous;
setup.functions.endpoint=@transferGAEndpoint;

% Auxdatas
setup.auxdata=auxdata;

% Boundaries
setup.bounds=bounds;

% Guess
setup.guess=guess;

% Solver
setup.nlp.solver='snopt';
setup.nlp.snoptoptions.maxiterations=1e3;

% Derivative
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';

% Mesh
setup.mesh.method = 'hp-LiuRao-Legendre';
setup.mesh.tolerance = 1e-9;
setup.mesh.maxiteration = 10;
setup.mesh.colpointsmax = 100;
setup.mesh.colpointsmin = 10;
for i=1:iphase
setup.mesh.phase(i).colpoints = 100*ones(1,10);
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
u=[u1];
unorm=sqrt(u(:,1).^2+u(:,2).^2+u(:,3).^2);
s1=output.result.solution.phase(1).state(:,1:6);
s=[s1];
for i=1:size(s1,1)
[r1(i,:),v1(i,:)]=mee2rv(s1(i,:),auxdata.mu_newUnit);
end
r=[r1];
x=r(:,1);
y=r(:,2);
z=r(:,3);

t1=output.result.solution.phase(1).time;
t=[t1];


subplot(1,2,1);
plot3(x,y,z,'LineWidth',1.5,'Color','black');hold on;
plot3(x(1),y(1),z(1),'c*','LineWidth',2);hold on;
text(x(1),y(1),z(1),'Departure');hold on;
plot3(x(end),y(end),z(end),'g*','LineWidth',2);hold on;
text(x(end),y(end),z(end),'Arrival');hold on;
plot3(0,0,0,'r*','LineWidth',2);hold on;
text(0,0,0,'Sun');hold on;

axis equal

subplot(1,2,2);
plot(t,unorm,'LineWidth',1.5,'Color','b');



















