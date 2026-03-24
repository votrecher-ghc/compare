% ============== parse_rinex_nav.m (同时支持GPS和北斗) ==============
function ephemeris = parse_rinex_nav_gps_bds(filepath)
% PARSE_RINEX_NAV - 解析RINEX 3.0x格式的导航文件，支持GPS和北斗(BDS)。
%
% 语法:  ephemeris = parse_rinex_nav(filepath)
%
% 输出:
%    ephemeris - 一个元胞数组，每个元胞的索引对应一个PRN号。
%                ephemeris{prn, sys_idx} 存储对应系统的星历。
%                sys_idx: 1 for GPS, 2 for BeiDou

% 为最多63颗卫星（BDS最多到C63）和2个系统初始化主元胞数组
% 第1列存GPS, 第2列存BeiDou
ephemeris = cell(63, 2);

% 打开文件进行读取
fid = fopen(filepath, 'r');
if fid == -1
    error('无法打开文件: %s', filepath);
end

% --- 读取文件头部分 ---
while true
    line = fgetl(fid);
    if contains(line, 'END OF HEADER')
        break;
    end
end

% --- 读取数据部分 ---
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || ~ischar(line) || length(line) < 80
        continue;
    end
    
    % 检查卫星系统标识符
    sys_char = upper(line(1));
    
    % ====================【核心修改：增加系统判断】====================
    if sys_char == 'G' || sys_char == 'C'
        
        % 为这个星历块创建一个新的结构体
        eph = struct();
        line = strrep(line, 'D', 'E'); % 统一科学计数法符号
        
        eph.PRN = str2double(line(2:3));
        eph.Toc.Year   = str2double(line(5:8));
        eph.Toc.Month  = str2double(line(10:11));
        eph.Toc.Day    = str2double(line(13:14));
        eph.Toc.Hour   = str2double(line(16:17));
        eph.Toc.Minute = str2double(line(19:20));
        eph.Toc.Second = str2double(line(22:23));
        
        switch sys_char
            case 'G' % --- 处理GPS (LNAV) 数据 ---
                eph.System = 'GPS';
                sys_idx = 1;
                
                % --- 第1行: af0, af1, af2 ---
                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));
                
                % 读取接下来的7行
                for i = 1:7
                    line = fgetl(fid);
                    line = strrep(line, 'D', 'E');
                    switch i
                        case 1; eph.IODE = str2double(line(5:23)); eph.Crs = str2double(line(24:42)); eph.Delta_n = str2double(line(43:61)); eph.M0 = str2double(line(62:80));
                        case 2; eph.Cuc = str2double(line(5:23)); eph.e = str2double(line(24:42)); eph.Cus = str2double(line(43:61)); eph.sqrtA = str2double(line(62:80));
                        case 3; eph.Toe = str2double(line(5:23)); eph.Cic = str2double(line(24:42)); eph.OMEGA0 = str2double(line(43:61)); eph.Cis = str2double(line(62:80));
                        case 4; eph.i0 = str2double(line(5:23)); eph.Crc = str2double(line(24:42)); eph.omega = str2double(line(43:61)); eph.OMEGA_DOT = str2double(line(62:80));
                        case 5; eph.IDOT = str2double(line(5:23)); eph.Codes_on_L2 = str2double(line(24:42)); eph.GPS_Week = str2double(line(43:61)); eph.L2_P_data_flag = str2double(line(62:80));
                        case 6; eph.SV_accuracy = str2double(line(5:23)); eph.SV_health = str2double(line(24:42)); eph.TGD = str2double(line(43:61)); eph.IODC = str2double(line(62:80));
                        case 7; eph.Transmission_time = str2double(line(5:23)); eph.Fit_interval = str2double(line(24:42));
                    end
                end

            case 'C' % --- 处理北斗 (BDS) 数据 ---
                eph.System = 'BeiDou';
                sys_idx = 2;

                % --- 第1行: A0, A1, A2 (北斗的钟差参数) ---
                eph.A0 = str2double(line(24:42));
                eph.A1 = str2double(line(43:61));
                eph.A2 = str2double(line(62:80));
                
                % 读取接下来的7行
                for i = 1:7
                    line = fgetl(fid);
                    line = strrep(line, 'D', 'E');
                    switch i
                        case 1; eph.IODE = str2double(line(5:23)); eph.Crs = str2double(line(24:42)); eph.Delta_n = str2double(line(43:61)); eph.M0 = str2double(line(62:80));
                        case 2; eph.Cuc = str2double(line(5:23)); eph.e = str2double(line(24:42)); eph.Cus = str2double(line(43:61)); eph.sqrtA = str2double(line(62:80));
                        case 3; eph.Toe = str2double(line(5:23)); eph.Cic = str2double(line(24:42)); eph.OMEGA0 = str2double(line(43:61)); eph.Cis = str2double(line(62:80));
                        case 4; eph.i0 = str2double(line(5:23)); eph.Crc = str2double(line(24:42)); eph.omega = str2double(line(43:61)); eph.OMEGA_DOT = str2double(line(62:80));
                        case 5; eph.IDOT = str2double(line(5:23)); eph.Spare1 = str2double(line(24:42)); eph.BDS_Week = str2double(line(43:61)); eph.Spare2 = str2double(line(62:80));
                        case 6; eph.SV_accuracy = str2double(line(5:23)); eph.SatH1 = str2double(line(24:42)); eph.TGD1 = str2double(line(43:61)); eph.TGD2 = str2double(line(62:80));
                        case 7; eph.Transmission_time = str2double(line(5:23)); eph.AODC = str2double(line(24:42));
                    end
                end
        end
        
        % 将新解析的星历结构体附加到元胞数组中
        if isempty(ephemeris{eph.PRN, sys_idx})
            ephemeris{eph.PRN, sys_idx} = eph;
        else
            ephemeris{eph.PRN, sys_idx}(end+1) = eph;
        end
    end
    % ======================================================================
end

% 关闭文件
fclose(fid);
fprintf('✅ 导航文件解析完成。\n');
end