% ============== calculate_average_sampling_rate.m ==============
function Fs = calculate_average_sampling_rate(obs_data)
% CALCULATE_AVERAGE_SAMPLING_RATE - 从解析后的观测数据中计算平均采样率。
%
% 语法: Fs = calculate_average_sampling_rate(obs_data)
%
% 输入:
%    obs_data - 由 parse_rinex_obs_... 函数解析出的结构体数组。
%
% 输出:
%    Fs       - 计算出的平均采样率 (Hz)。

fprintf('--> 正在计算平均采样率...\n');

% 检查输入数据点是否足够
if length(obs_data) < 2
    error('观测数据点不足 (少于2个)，无法计算采样率。');
end

% 从结构体数组中提取所有的时间戳，构成一个datetime向量
time_vector = [obs_data.time];

% 计算相邻时间戳之间的平均时间差（单位：秒）
mean_dt = mean(seconds(diff(time_vector)));

% 采样率是平均时间间隔的倒数
Fs = 1 / mean_dt;

fprintf('    计算出的平均采样率为: %.2f Hz\n', Fs);
fprintf('✅ 平均采样率计算完成。\n\n');

end