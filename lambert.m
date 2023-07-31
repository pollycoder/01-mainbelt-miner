



function [v1, v2] = lambert(r1, r2, dt, mu, is_prograde)
    % r1: 初始位置向量 [x; y; z]
    % r2: 目标位置向量 [x; y; z]
    % dt: 转移时间 (s)
    % mu: 万有引力常数 (m^3/s^2)
    % is_prograde: 指示是否使用顺行法则 (true/false)

    % 求解Lambert问题
    [v1, v2] = lambert(r1, r2, dt, mu, is_prograde, 0, 0);

    % 调整速度向量的方向
    if is_prograde
        v1 = v1;
        v2 = v2;
    else
        v1 = -v1;
        v2 = -v2;
    end
end
