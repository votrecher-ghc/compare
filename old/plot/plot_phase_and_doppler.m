% ============== plot_phase_and_doppler.m ==============
function plot_phase_and_doppler(obs_data, satellite_id, obs_code_suffix)
% PLOT_PHASE_AND_DOPPLER - 对比绘制指定卫星的原始载波相位和多普勒观测值。
%
% 语法: plot_phase_and_doppler(obs_data, satellite_id, obs_code_suffix)
%
% 输入:
%    obs_data          - 由 parse_rinex_obs_advanced 解析出的观测数据。
%    satellite_id      - 目标卫星ID (例如, 'G10')。
%    obs_code_suffix   - 观测码的后缀 (例如, '1C', '2W' 等)。
%                        函数会自动匹配对应的L码和D码。

% --- 1. 参数检查与初始化 ---
if nargin < 3
    error('输入参数不足。请提供观测数据、卫星ID和观测码后缀。');
end

% 根据输入后缀，构建完整的载波相位和多普勒观测码
full_phase_code = ['L', obs_code_suffix];
full_doppler_code = ['D', obs_code_suffix];

num_epochs = length(obs_data);
% 初始化用于存储结果的数组
time_vector = NaT(0,0); 
phase_rad_vector = [];  % 存储相位 (弧度)
doppler_hz_vector = []; % 新增：用于存储多普勒 (Hz)

fprintf('--> 开始提取卫星 %s 的 %s 和 %s 原始数据...\n', satellite_id, full_phase_code, full_doppler_code);

% --- 2. 遍历所有历元，提取并转换数据 ---
for i = 1:num_epochs
    % 检查当前历元是否观测到了目标卫星
    if isfield(obs_data(i).data, satellite_id)
        sat_data = obs_data(i).data.(satellite_id);
        
        % 检查该卫星是否同时包含指定的相位和多普勒观测值
        if isfield(sat_data, 'carrier_phase') && isfield(sat_data.carrier_phase, full_phase_code) && ...
           isfield(sat_data, 'doppler') && isfield(sat_data.doppler, full_doppler_code)
            
            phase_cycles = sat_data.carrier_phase.(full_phase_code);
            doppler_hz = sat_data.doppler.(full_doppler_code);
            
            % 只有当两个值都有效时，才进行存储
            if ~isnan(phase_cycles) && ~isnan(doppler_hz)
                % 存储时间
                time_vector(end+1,1) = obs_data(i).time;
                % 存储多普勒值 (单位: Hz)
                doppler_hz_vector(end+1,1) = doppler_hz;
                % 将相位从“周”转换为“弧度”并存储
                phase_rad_vector(end+1,1) = phase_cycles * (2 * pi);
            end
        end
    end
end

fprintf('--> 数据提取完成！共找到 %d 个有效的成对数据点。\n\n', length(time_vector));

% --- 3. 数据检查与绘图 ---
if isempty(time_vector)
    fprintf('警告: 在所有历元中，均未找到卫星 %s 的 %s 和 %s 成对观测数据，无法绘图。\n', satellite_id, full_phase_code, full_doppler_code);
    return;
end

% 创建一个新的图形窗口
figure('Name', sprintf('卫星 %s 的原始相位与多普勒对比', satellite_id), 'NumberTitle', 'off');

% --- 子图1: 原始的多普勒频移 ---
ax1 = subplot(2, 1, 1);
plot(time_vector, doppler_hz_vector, 'b.-');
title(sprintf('卫星 %s 原始多普勒频移 (%s)', satellite_id, full_doppler_code));
xlabel('时间');
ylabel('多普勒 (Hz)');
grid on;
datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');

% --- 子图2: 转换为“弧度”并展开后的相位 ---
ax2 = subplot(2, 1, 2);
plot(time_vector, unwrap(phase_rad_vector), 'r.-');
title(sprintf('卫星 %s 原始累积载波相位 (%s)', satellite_id, full_phase_code));
xlabel('时间');
ylabel('累积相位 (弧度)');
grid on;
datetick('x', 'HH:MM:SS', 'keepticks', 'keeplimits');

% 联动两个子图的X轴，方便同步缩放和查看
linkaxes([ax1, ax2], 'x');

fprintf('✅ 图形绘制完成！\n');

end