function E = f2E(f, e)
    if e < 0.0
        E = f;
        return;
    end
    
    E = 0.0;
    
    if e >= 0.0 && e < 1.0 % 圆和椭圆轨道
        f0 = mod(f, 2*pi);
        if f0 > pi
            f0 = f0 - 2*pi;
        end
        if f0 < -pi
            f0 = f0 + 2*pi;
        end
        E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(0.5 * f0));
        E = E + f - f0;
    elseif e > 1.0 % 双曲线轨道
        if (f > pi - acos(1.0 / e)) || (f < -pi + acos(1.0 / e))
            % disp('不可能达到的双曲轨道.');
            E = f;
            return;
        else
            E = 2.0 * atanh(sqrt((e - 1.0) / (1.0 + e)) * tan(0.5 * f));
        end
    else % 抛物线轨道
        E = f;
        % disp('抛物线轨道没有定义偏近点角.在此将其值置为f.');
    end
end
