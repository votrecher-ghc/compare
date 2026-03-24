% --- 使用示例 ---
clear;
clc;
close all;

%% 1. 解析导航文件
% 请将 'COM5___460800_250614_122524.24n' 替换为您的实际文件路径
nav_filepath = 'COM5___460800_250616_032142.nav'; 
ephemeris_data = parse_rinex_nav(nav_filepath);

%% 2. 定义观测时间和目标卫星
obs_time_gps_seconds = 568810; % 示例: GPS周内秒
obs_week = 2370;               % 示例: GPS周
prn_to_find = 9;              % 我们要计算G10卫星的位置

%% 3. 寻找给定时间的最佳星历
eph_sets = ephemeris_data{prn_to_find};
if isempty(eph_sets)
    error('在导航文件中未找到PRN %d的任何星历数据。', prn_to_find);
end

% 寻找 |t_obs - Toe| 最小的星历
min_diff = inf;
best_eph_index = -1;

for i = 1:length(eph_sets)
    % 注意: 完整的实现需要将观测时间转换为与Toe一致的时间尺度，
    % 并考虑周翻转。为简化起见，此处假设它们在同一周内。
    time_diff = abs(obs_time_gps_seconds - eph_sets(i).Toe);
    if time_diff < min_diff
        min_diff = time_diff;
        best_eph_index = i;
    end
end

% 检查是否找到了有效的星历 (通常在2-4小时窗口内)
% 修正后的判断语句
if best_eph_index == -1 || min_diff > 4*3600 
    error('在指定时间未找到适用于PRN %d的星历。', prn_to_find);
end

%% 4. 访问并显示选定的星历参数
best_eph = eph_sets(best_eph_index);
eccentricity = best_eph.e;
toe = best_eph.Toe;

fprintf('为G%d找到最佳星历，其Toe = %.1f s\n', prn_to_find, toe);
fprintf('轨道偏心率 (e) = %e\n', eccentricity);