%%%%%%%%%%%%%%%%%%%%%%%%%
% SOI
% Input: 
%   r:3x1 vector
%   v:3x1 vector
%%%%%%%%%%%%%%%%%%%%%%%%%
function [rRel_out,vRel_out,time]=SOI(mu,rSOI,rRel_in,vRel_in)
cosine_in=abs((rRel_in'*vRel_in)/(norm(rRel_in,2)*norm(vRel_in,2)));
sine_in=sqrt(1-cosine_in.^2);
d=rSOI*sine_in;

theta=2.*asin(mu./(norm(vRel_in,2).^2.*d));
rotate_v=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
vRel_out=rotate_v*vRel_in;

alpha=acos(d./rSOI)+theta./2;
alpha=2*alpha;
rotate_r=[cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
rRel_out=rotate_r*rRel_in;

coe1=abs(rv2coe(rRel_in,vRel_in,mu));
coe2=abs(rv2coe(rRel_out,vRel_out,mu));
time=f0ft2dt(coe1(6),coe2(6),coe1(1),coe1(2),mu);
end