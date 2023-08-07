function [r1,v1] = rv02rvf(r0,v0, dt, mu)
    coe = rv2coe(r0,v0, mu);
    
    coe(6) = f0dt2ft(coe(6), dt, coe(1), coe(2), mu);
    coe(6) = mod(coe(6), 2*pi);
    [r1,v1] = coe2rv(coe, mu);

end
