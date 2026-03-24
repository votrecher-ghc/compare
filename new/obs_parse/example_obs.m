clear;
clc;
close all;


% 解析观测文件
obs_data = parse_rinex_obs_advanced('COM5___460800_250614_122524.obs');
disp(obs_data(1).time)  
% 解析导航文件
% nav_data = parse_rinex_nav('COM5___460800_250614_122524.24n');

%%
% 画出出某个卫星的解析数值
plot_observations_per_satellite(obs_data, 'G10', '1C');  % GPS L1C
%%
% 列出观测的卫星
list_observed_satellites(obs_data);

%%
%打印出每个系统的观测码
print_observation_types_by_system(obs_data);

%%
[time, bias] = estimate_receiver_clock_bias(obs_data, 'G', 'C1C');

%%
[t,phi,s] = reconstruct_phase_with_code(obs_data,'G10','1C');   % GPS L1C

%%
plot_phase_and_signal(obs_data, 'G10', '1C');