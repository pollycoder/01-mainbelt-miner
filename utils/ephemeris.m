function [pos, vel] = ephemeris(planet) 

    % 设置JPL DE421星历文件路径
    % 根据实际情况修改路径
    ephemeris_file_path = './ephemeris/de421.bsp';

    % 加载星历文
    cspice_furnsh('./mice/data/cook_01.tls');
    cspice_furnsh(ephemeris_file_path);

    % 定义时间范围（示例：从2023年1月1日到2023年12月31日）
    start_time = '2023 JAN 01 00:00:00.000';
    end_time = '2030 DEC 31 23:59:59.999';

    % 设定时间间隔（单位：秒，示例设置为1天）
    time_step = cspice_spd;

    % 生成时间序列
    times = cspice_str2et({start_time, end_time});
    et_time_series = (times(1):time_step:times(2));

    % 初始化位置和速度数组
    num_steps = numel(et_time_series);
    pos = zeros(3, num_steps);
    vel = zeros(3, num_steps);

    % 逐步获取位置和速度数据
    for i = 1:num_steps
        [state, ~] = cspice_spkezr(planet, et_time_series(i), 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER');
       
        pos(:, i) = state(1:3);
        vel(:, i) = state(4:6);
    end

    % 卸载星历文件
    cspice_kclear;
end
