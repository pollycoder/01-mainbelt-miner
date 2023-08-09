function coe = rv2coe(r, v, mu)
    % 输入：位置向量r，速度向量v，引力常数mu
    % 输出：计算得到的轨道根数coe
    
    coe = zeros(6, 1);
    
    r_norm = norm(r);
    v_norm = norm(v);
    
    h_vec = cross(r, v);
    h = norm(h_vec);
    
    % 计算半通径
    a = 1 / (2 / r_norm - v_norm^2 / mu);
    coe(1) = a;
    
    % 计算偏心率向量
    e_vec = cross(v, h_vec) / mu - r / r_norm;
    e = norm(e_vec);
    coe(2) = e;
    
    % 计算倾角
    inc = acos(h_vec(3) / h);
    coe(3) = inc;
    
    % 计算升交点赤经
    N = cross([0; 0; 1], h_vec);
    N_norm = norm(N);
    RAAN = acos(N(1) / N_norm);
    if N(2) < 0
        RAAN = 2 * pi - RAAN;
    end
    coe(4) = RAAN;
    
    % 计算近星点幅角、真近点角或偏近点角
    if e < 1.0 - 1e-10  % 椭圆轨道
        E = atan2(dot(e_vec, cross(h_vec, r)) / (e * h), dot(r, e_vec) / mu - 1 / r_norm);
        coe(5) = E;
    elseif e > 1.0 + 1e-10  % 双曲线轨道
        F = acosh((dot(e_vec, r) / (e * r_norm) - 1) / e);
        coe(5) = F;
    else  % 抛物线轨道
        B = cross(h_vec, r) / h;
        coe(5) = atan2(real(dot(B, v)) / sqrt(mu), dot(r, v) / sqrt(mu));
    end
    
    % 计算真近点角或偏近点角对应的时间参数
    if e < 1.0  % 椭圆轨道
        M = E - e * sin(E);
    elseif e > 1.0  % 双曲线轨道
        M = e * sinh(F) - F;
    else  % 抛物线轨道
        M = coe(5);
    end
    coe(6) = M;
end
