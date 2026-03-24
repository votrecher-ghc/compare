function plot_observations_per_satellite(data, sat_id, obs_band)
% 绘制某颗卫星的 C/L/D/S 四种观测值时间序列图（假设 30Hz 采样率）
%
% 示例：
% plot_observations_per_satellite(data, 'G10', '1C')

    % --- 1. 构造观测码 ---
    codeC = ['C' obs_band];
    codeL = ['L' obs_band];
    codeD = ['D' obs_band];
    codeS = ['S' obs_band];
    
    % --- 2. 初始化时间 & 观测值容器 ---
    t_all = datetime.empty(0,1);
    C = []; L = []; D = []; S = [];
    
    % --- 3. 遍历历元提取数据 ---
    for i = 1:numel(data)
        if isfield(data(i).data, sat_id)
            obs = data(i).data.(sat_id);
            t_all(end+1,1) = data(i).time;
            C(end+1,1) = get_obs_val(obs, 'pseudorange', codeC);
            L(end+1,1) = get_obs_val(obs, 'carrier_phase', codeL);
            D(end+1,1) = get_obs_val(obs, 'doppler', codeD);
            S(end+1,1) = get_obs_val(obs, 'snr', codeS);
        end
    end
    
    % 如果没有找到该卫星的数据，则直接返回
    if isempty(t_all)
        disp(['警告: 未在数据中找到卫星 ' sat_id ' 的观测值。']);
        return;
    end

    % --- 4. 构造均匀时间轴并插值 ---
    % 仅当有多个数据点时才进行插值
    if numel(t_all) > 1
        t_uniform = t_all(1):seconds(1/30):t_all(end);
        
        % 使用 'linear' 方法并用 NaN 填充外插值
        C = interp1(t_all, C, t_uniform, 'linear', NaN);
        L = interp1(t_all, L, t_uniform, 'linear', NaN);
        % >> 修正点: 此处添加了 'linear' 方法 <<
        D = interp1(t_all, D, t_uniform, 'linear', NaN);
        S = interp1(t_all, S, t_uniform, 'linear', NaN);
    else
        % 如果只有一个点，则无需插值
        t_uniform = t_all;
    end
    
    % --- 5. 绘图 (使用 Subplot) ---
    fig = figure('Name',['观测值 - ' sat_id ' (' obs_band ')'], 'Position', [100 100 900 700]);
    
    ax1 = subplot(4,1,1);
    plot(t_uniform, C, 'b'); ylabel(codeC); grid on;
    title(['伪距 ' codeC]);
    
    ax2 = subplot(4,1,2);
    plot(t_uniform, L, 'g'); ylabel(codeL); grid on;
    title(['载波相位 ' codeL]);
    
    ax3 = subplot(4,1,3);
    plot(t_uniform, D, 'r'); ylabel(codeD); grid on;
    title(['多普勒 ' codeD]);
    
    ax4 = subplot(4,1,4);
    plot(t_uniform, S, 'm'); ylabel(codeS); grid on;
    xlabel('时间'); title(['信噪比 ' codeS]);
    
    % >> 改进点: 联动所有子图的X轴，方便缩放和查看 <<
    linkaxes([ax1, ax2, ax3, ax4], 'x');
end

% 辅助函数，保持不变
function val = get_obs_val(obs, field, code)
    if isfield(obs, field) && isfield(obs.(field), code)
        val = obs.(field).(code);
    else
        val = NaN;
    end
end