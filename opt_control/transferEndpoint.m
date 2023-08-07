%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint function - all new unit
% No input.
% Object:
%   The integrant
%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=transferEndpoint(input)
% Earth to Mars
t01=input.phase(1).initialtime;
tf1=input.phase(1).finaltime;
x01=input.phase(1).initialstate;
xf1=input.phase(1).finalstate;

% Mars to Main Belt
t02=input.phase(2).initialtime;
tf2=input.phase(2).finaltime;
x02=input.phase(2).initialstate;
xf2=input.phase(2).finalstate;

% GA-Mars-120°
% Mars-coordinate system
rRelative1_newUnit=x01(1:3)-input.auxdata.MarsState_newUnit(1:3);
vRelative1_newUnit=x01(4:6)-input.auxdata.MarsState_newUnit(4:6);

thetav=pi/3;
thetar=2*pi/3;
rotate_v=[cos(thetav),-sin(thetav),0;sin(thetav),cos(thetav),0;0,0,1];
rotate_r=[cos(thetar),-sin(thetar),0;sin(thetar),cos(thetar),0;0,0,1];
vRelative2_newUnit=(rotate_v*vRelative1_newUnit')';
rRelative2_newUnit=(rotate_r*rRelative1_newUnit')';

% Transfer back to heliocentric system
v0GA_newUnit=vRelative2_newUnit+input.auxdata.MarsState_newUnit(4:6);
r0GA_newUnit=rRelative2_newUnit+input.auxdata.MarsState_newUnit(1:3);
x0GA_newUnit=[r0GA_newUnit,v0GA_newUnit,x01(7)];

output.objective=-xf2(7);
output.eventgroup.event=[x02-x0GA_newUnit,t02-tf1];
end