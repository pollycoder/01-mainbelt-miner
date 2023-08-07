function ft = f0dt2ft(f0, dt, a, e, mu, MaxIter, epsilon)
    if nargin==5
        MaxIter=100;
        epsilon=1e-6;
    end
    if mu <= 0.0 || MaxIter < 1 || a <= 0.0 || e < 0.0
        ft = f0;
        return;
    end
    
    ft = 0.0;
    
    if (e >= 0.0 && e < 1.0) || (e > 1.0) % 圆、椭圆、双曲轨道
        E = f2E(f0, e);
        M = E2M(E, e);
        M = M + sqrt(mu / (a^3)) * dt;
        E = M2E(M, e, MaxIter, epsilon);
        ft = E2f(E, e);
    else % 抛物线轨道
        if f0 < -pi || f0 > pi
            % disp('对于抛物线轨道，初始真近点角应在-180至180度之间.');
           
            ft = f0;
            return;
        elseif f0 > pi || f0 < -pi
            ft = f0;
        else
            B = 0.75 * sqrt(2.0 * mu / (a^3)) * dt + 0.5 * tan(0.5 * f0) * ((tan(0.5 * f0))^2 + 3.0);
            B1B = B + sqrt(1.0 + B^2);
            tanv = 0.0;
            if abs(dt) < 2*pi * sqrt((a^3) / mu) / 1000.0 % 推进时间为小量的情况
                A = (B1B)^(2/3);
                tanv = 2.0 * A * B / (1.0 + (1.0 + A) * A);
            else % 不是小量的情况
                temp = (B1B)^(1/3);
                tanv = temp - 1.0 / temp;
            end
            ft = 2.0 * atan(tanv);
        end
    end
end
