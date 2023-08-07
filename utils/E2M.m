function M = E2M(E, e)
    if e < 0.0
        M = E;
        return;
    end
    
    if e >= 0.0 && e < 1.0 % 圆和椭圆轨道
        E0 = mod(E, 2*pi);
        M = E0 - e * sin(E0);
        M = M + E - E0;
    elseif e > 1.0 % 双曲线轨道
        M = e * sinh(E) - E;
    else % 抛物线轨道
        M = E;
        % disp('抛物线轨道没有定义偏近点角.在此将平近点角值置为E.');
    end
end
