function [p,f,g,h,k,L] = coe2mee(a,e,i,Ome,ome,theta)
p = a*(1-e^2);
f = e*cos(ome+Ome);
g = e*sin(ome+Ome);
h = tan(i/2)*cos(Ome);
k = tan(i/2)*sin(Ome);
L = Ome+ome+theta;
end