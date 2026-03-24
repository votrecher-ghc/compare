function unique_sats = list_observed_satellites(data)
% 列出 RINEX 数据中所有观测到的唯一卫星编号（如 G01, C29, etc.）
% 以每行约8个卫星的格式输出

    all_sats = {};

    for i = 1:length(data)
        if isfield(data(i), 'data')
            sats = fieldnames(data(i).data);
            all_sats = [all_sats; sats(:)];
        end
    end

    unique_sats = unique(all_sats);
    num_sats = length(unique_sats);
    
    % 输出结果
    fprintf('✅ 共观测到 %d 颗不同卫星：\n', num_sats);
    
    % 设置每行显示的卫星数量
    sats_per_line = 8;
    
    % 计算需要显示的行数
    num_lines = ceil(num_sats / sats_per_line);
    
    % 格式化输出卫星列表
    for i = 1:num_lines
        start_idx = (i-1)*sats_per_line + 1;
        end_idx = min(i*sats_per_line, num_sats);
        
        % 创建当前行的格式字符串
        line_format = repmat('%-6s', 1, end_idx - start_idx + 1);
        
        % 输出当前行的卫星
        fprintf([line_format '\n'], unique_sats{start_idx:end_idx});
    end
end