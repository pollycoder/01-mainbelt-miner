


function [global_best, global_best_position] = pso(obj_func, num_particles, num_iterations, lower_bound, upper_bound)
    % obj_func: 优化目标函数的句柄，形如 @(x) func(x)，其中x是输入参数，func是目标函数
    % num_particles: 粒子数量
    % num_iterations: 迭代次数
    % lower_bound: 参数的下界
    % upper_bound: 参数的上界

    % 初始化粒子群
    num_dimensions = numel(lower_bound);
    particles_position = rand(num_particles, num_dimensions) .* (upper_bound - lower_bound) + lower_bound;
    particles_velocity = zeros(num_particles, num_dimensions);
    particles_best_position = particles_position;
    particles_best_value = arrayfun(obj_func, particles_position);

    % 初始化全局最优
    [global_best_value, global_best_index] = min(particles_best_value);
    global_best = global_best_value;
    global_best_position = particles_best_position(global_best_index, :);

    % PSO主循环
    for iter = 1:num_iterations
        for i = 1:num_particles
            % 更新速度和位置
            inertia_weight = 0.8; % 惯性权重
            cognitive_weight = 2.0; % 认知权重
            social_weight = 2.0; % 社会权重

            r1 = rand(1, num_dimensions);
            r2 = rand(1, num_dimensions);

            particles_velocity(i, :) = inertia_weight * particles_velocity(i, :) ...
                + cognitive_weight * r1 .* (particles_best_position(i, :) - particles_position(i, :)) ...
                + social_weight * r2 .* (global_best_position - particles_position(i, :));

            % 限制速度范围
            max_velocity = (upper_bound - lower_bound) * 0.2;
            particles_velocity(i, :) = max(particles_velocity(i, :), -max_velocity);
            particles_velocity(i, :) = min(particles_velocity(i, :), max_velocity);

            % 更新位置
            particles_position(i, :) = particles_position(i, :) + particles_velocity(i, :);

            % 限制位置范围
            particles_position(i, :) = max(particles_position(i, :), lower_bound);
            particles_position(i, :) = min(particles_position(i, :), upper_bound);

            % 更新个体最优
            new_value = obj_func(particles_position(i, :));
            if new_value < particles_best_value(i)
                particles_best_value(i) = new_value;
                particles_best_position(i, :) = particles_position(i, :);

                % 更新全局最优
                if new_value < global_best_value
                    global_best_value = new_value;
                    global_best = global_best_value;
                    global_best_position = particles_best_position(i, :);
                end
            end
        end

        fprintf('Iteration %d: Best Value = %.4f\n', iter, global_best_value);
    end
end
