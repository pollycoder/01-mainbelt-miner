function E = M2E(M, e, MaxIter, epsilon)
    if nargin==2
        MaxIter=100;
        epsilon=1e-6;
    end
    Delta3 = 0.0;
    N = 0;
    
    if e >= 0.0 && e < 1.0 % 圆和椭圆轨道
        RM = mod(M, 2*pi);
        if RM < 0.0
            RM = RM + 2*pi;
        end
        sinRM = sin(RM);
        E = RM + 0.85 * e * sign(sinRM);
        Delta3 = 1.0;
        while abs(Delta3) >= epsilon && N < MaxIter
            Minus = E - e * sin(E) - RM;
            DeMinus = 1.0 - e * cos(E);
            DeDeMinus = e * sin(E);
            DeDeDeMinus = e * cos(E);
            Delta1 = -Minus / DeMinus;
            Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
            Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
            E = E + Delta3;
            N = N + 1;
        end
        E = E + M - RM;
    elseif e > 1.0 % 双曲线轨道
        E = asinh(M / e);
        Delta3 = 1.0;
        while abs(Delta3) >= epsilon && N < MaxIter
            Minus = e * sinh(E) - E - M;
            DeMinus = e * cosh(E) - 1.0;
            DeDeMinus = e * sinh(E);
            DeDeDeMinus = e * cosh(E);
            Delta1 = -Minus / DeMinus;
            Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
            Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
            E = E + Delta3;
            N = N + 1;
        end
    else % 抛物线轨道
        E = M;
        % disp('抛物线轨道没有定义偏近点角.在此将其值置为M.');
    end
    
    if (((e >= 0.0 && e < 1.0) || (e > 1.0)) && abs(Delta3) >= 5.0 * epsilon && N >= MaxIter)
        % disp('迭代不收敛,请降低精度epsilon或增加迭代次数限制.');
        E = M;
        return;
    end
end
