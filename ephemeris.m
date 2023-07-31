function [earth_pos, earth_vel, mars_pos, mars_vel] = ephemeris()

    % 设置JPL DE421星历文件路径
    % 根据实际情况修改路径
    ephemeris_file_path = 'path/to/DE421.bsp';

    % 加载星历文件
    spice_toolbox = '../mice/';  % 请替换为你的SPICE工具箱路径
    addpath(spice_toolbox);
    cspice_furnsh(ephemeris_file_path);

    % 定义时间范围（示例：从2023年1月1日到2023年12月31日）
    start_time = '2023 JAN 01 00:00:00.000';
    end_time = '2023 DEC 31 23:59:59.999';

    % 设定时间间隔（单位：秒，示例设置为1天）
    time_step = cspice_spd;

    % 生成时间序列
    times = cspice_str2et({start_time, end_time});
    et_time_series = (times(1):time_step:times(2));

    % 初始化位置和速度数组
    num_steps = numel(et_time_series);
    earth_pos = zeros(3, num_steps);
    earth_vel = zeros(3, num_steps);
    mars_pos = zeros(3, num_steps);
    mars_vel = zeros(3, num_steps);

    % 逐步获取位置和速度数据
    for i = 1:num_steps
        [earth_state, ~] = cspice_spkezr('EARTH', et_time_series(i), 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER');
        [mars_state, ~] = cspice_spkezr('MARS', et_time_series(i), 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER');

        earth_pos(:, i) = earth_state(1:3);
        earth_vel(:, i) = earth_state(4:6);
        mars_pos(:, i) = mars_state(1:3);
        mars_vel(:, i) = mars_state(4:6);
    end

    % 卸载星历文件
    cspice_kclear;
end
