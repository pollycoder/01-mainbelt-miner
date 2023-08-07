function [a,e,i,Ome,ome,theta] = mee2coe(p,f,g,h,k,L)
%半长轴、离心率、轨道倾角、升交点经度、近心点辐角、真近点角;
a = p/(1-f*f-g*g);
e = sqrt(f*f+g*g);
i = atan2(2*sqrt(h*h+k*k),1-h*h-k*k);
ome = atan2(g*h-f*k,f*h+g*k);
Ome = atan2(k,h);
theta = L-Ome-ome;
end