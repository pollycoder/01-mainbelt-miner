function [r,v]=mee2rv(mee,mu)
p=mee(1);
f=mee(2);
g=mee(3);
h=mee(4);
k=mee(5);
L=mee(6);
q=1+f.*cos(L)+g.*sin(L);
r0=p./q;
alpha2=h.^2-k.^2;
s2=1+h.^2+k.^2;
r=r0./s2.*[cos(L)+alpha2.*cos(L)+2.*h.*k.*sin(L); 
          sin(L)-alpha2.*sin(L)+2.*h.*k.*cos(L); 
          2.*(h.*sin(L)-k.*sin(L))];
v=-1/s2.*sqrt(mu/p).*[sin(L)+alpha2.*sin(L)-2.*h.*k.*cos(L)+g-2.*f.*k.*h+alpha2.*g; 
                   -cos(L)+alpha2.*cos(L)-2.*h.*k.*sin(L)-f+2.*g.*k.*h+alpha2.*f;
                   -(h.*cos(L)+k.*sin(L)+f.*h+g.*k)];
end