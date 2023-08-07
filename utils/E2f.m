function f = E2f(E, e)
    if e < 0.0
        f = E;
        return;
    end
    
    if e >= 0.0 && e < 1.0 % 圆和椭圆轨道
        E0 = mod(E, 2*pi);
        if E0 > pi
            E0 = E0 - 2*pi;
        end
        if E0 < -pi
            E0 = E0 + 2*pi;
        end
        f = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(0.5 * E0));
        f = f + E - E0;
    elseif e > 1.0 % 双曲线轨道
        f = 2.0 * atan(sqrt((e + 1.0) / (e - 1.0)) * tanh(0.5 * E));
    else % 抛物线轨道
        f = E;
        % disp('抛物线轨道没有定义偏近点角.在此将真近点角的值置为E.');
    end
end
