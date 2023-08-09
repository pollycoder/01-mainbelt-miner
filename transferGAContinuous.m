%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous function - all new unit
%%%%%%%%%%%%%%%%%%%%%%%%%

function output=transferGAContinuous(input)
Isp_newUnit=input.auxdata.Isp_newUnit;                         
g0_newUnit=input.auxdata.g0_newUnit;
Tmax_newUnit=input.auxdata.Tmax_newUnit;
muSun=input.auxdata.mu_newUnit;
muMars=input.auxdata.muMars_newUnit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first phase - Earth--->Mars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iphase=1;
x_newUnit=input.phase(iphase).state;
u_newUnit=input.phase(iphase).control;

% The elements
p=x_newUnit(:,1);
f=x_newUnit(:,2);
g=x_newUnit(:,3);
h=x_newUnit(:,4);
k=x_newUnit(:,5);
L=x_newUnit(:,6);
m=x_newUnit(:,7);

% Points
num_points=length(p);
dynamics=zeros(7,num_points);
for i=1:num_points
    W=1+f(i)*cos(L(i))+g(i)*sin(L(i));
    Z=h(i)*sin(L(i))-k(i)*cos(L(i));
    C=1+h(i)^2+k(i)^2;
    A=(sqrt(p(i)/muSun)/W)*[0;0;0;0;0;W^3*muSun/(p(i)^2)];
    B=sqrt(p(i)/muSun)/W* ...
        [0,            2*p(i),                 0;
         W*sin(L(i)),  (W+1)*cos(L(i))+f(i),  -Z*g(i);
        -W*cos(L(i)),  (W+1)*sin(L(i))+g(i),   Z*f(i); 
         0,            0,                      1/2*C*cos(L(i));
         0,            0,                      1/2*C*sin(L(i));
         0,            0,                      Z];

    acc_eng=Tmax_newUnit*u_newUnit(i,:)/m(i);
    dynamics(1:6,i)=B*acc_eng'+A;

    unorm(i)=sqrt(u_newUnit(i,1)^2+u_newUnit(i,2)^2+u_newUnit(i,3)^2);
    dynamics(7,i)=-Tmax_newUnit*unorm(i)/(Isp_newUnit*g0_newUnit);
end
output(iphase).dynamics=dynamics';
output(iphase).path=unorm;

end
