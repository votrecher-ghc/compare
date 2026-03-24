function print_observation_types_by_system_fast(data)
% 更高效地统计各导航系统的观测码类型（如 C1C, L1C 等）
% 输入：data —— 由 parse_rinex_obs_geopp 返回的结构数组

    system_map = struct();  % 每个系统一个 map 存储唯一观测码
    label_map = struct('G','GPS','C','BDS','E','Galileo','R','GLONASS',...
                       'J','QZSS','I','IRNSS','S','SBAS');

    for i = 1:numel(data)
        sats = fieldnames(data(i).data);
        for s = 1:numel(sats)
            sat_id = sats{s};
            sys = sat_id(1);  % 'G', 'C', etc.
            obs = data(i).data.(sat_id);

            if ~isfield(system_map, sys)
                system_map.(sys) = containers.Map();
            end

            fields = {'pseudorange','carrier_phase','doppler','snr'};
            for f = 1:length(fields)
                fld = fields{f};
                if isfield(obs, fld)
                    codes = fieldnames(obs.(fld));
                    for c = 1:numel(codes)
                        code = codes{c};
                        system_map.(sys)(code) = true;  % 加入 map，自动去重
                    end
                end
            end
        end
    end

    % ---------- 打印结果 ----------
    fprintf('\n各导航系统观测码统计：\n');
    sys_keys = fieldnames(system_map);
    for i = 1:numel(sys_keys)
        key = sys_keys{i};
        code_list = keys(system_map.(key));
        code_list = sort(code_list);  % 排序更美观
        sys_name = key;
        if isfield(label_map, key)
            sys_name = label_map.(key);
        end
        fprintf('%-8s 共 %2d 种观测码：%s\n', sys_name, numel(code_list), strjoin(code_list, ', '));
    end
end
