function dt = f0ft2dt(f0, ft, a, e, mu)
    if mu <= 0.0 || a <= 0.0 || e < 0.0
        dt = 0.0;
        return;
    end
    
    dt = 0.0;
    if e >= 1.0
        maxangle = pi - acos(1.0 / e);
        if (min(f0, ft) < -maxangle) || (max(f0, ft) > maxangle)
            % 不可能达到的双曲或抛物线轨道
            dt = 0.0;
            return;
        elseif (f0 < -maxangle || f0 > maxangle || ft < -maxangle || ft > maxangle)
            % 所需时间很难计算准确，可能为无穷，在此置为一个大数
            dt = 1.0e308;
            return;
        end
    end

    omega = sqrt(mu / (a * a * a));
    delta = 0.0;
    if (e >= 0.0 && e < 1.0) || (e > 1.0)
        E_f0 = f2E(f0, e);
        M0 = E2M(E_f0, e);
        E_ft = f2E(ft, e);
        Mt = E2M(E_ft, e);
        delta = Mt - M0;
    else
        B1 = tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
        B2 = tan(0.5 * ft) * ((tan(0.5 * ft)) * (tan(0.5 * ft)) + 3.0);
        delta = sqrt(2.0) / 3.0 * (B2 - B1);
    end
    dt = delta / omega;
end