function [time_vec, value_vec] = extract_satellite_observation(data, target_sat, obs_type)
    % 提取指定卫星、指定观测码的时间和值
    
    time_vec = datetime.empty(0, 1);
    value_vec = [];

    for i = 1:length(data)
        if isfield(data(i).data, target_sat)
            sat_data = data(i).data.(target_sat);
            val = NaN;

            % 分类查找观测码
            if isfield(sat_data, 'pseudorange') && isfield(sat_data.pseudorange, obs_type)
                val = sat_data.pseudorange.(obs_type);
            elseif isfield(sat_data, 'carrier_phase') && isfield(sat_data.carrier_phase, obs_type)
                val = sat_data.carrier_phase.(obs_type);
            elseif isfield(sat_data, 'doppler') && isfield(sat_data.doppler, obs_type)
                val = sat_data.doppler.(obs_type);
            elseif isfield(sat_data, 'snr') && isfield(sat_data.snr, obs_type)
                val = sat_data.snr.(obs_type);
            end

            % 收集非缺失值
            if ~isnan(val)
                time_vec(end+1, 1)  = data(i).time;
                value_vec(end+1, 1) = val;
            end
        end
    end
end
