% ============== parse_rinex_nav_multi_gnss.m (G=1, C=2, R=3, E=4, J=5) ==============
function ephemeris = parse_rinex_nav_multi_gnss(filepath)
% PARSE_RINEX_NAV_MULTI_GNSS Parse RINEX 3.xx navigation file.
%
% Output layout (legacy compatible):
%   ephemeris{PRN,1} -> GPS (G)
%   ephemeris{PRN,2} -> BeiDou (C)
%   ephemeris{PRN,3} -> GLONASS (R)
%   ephemeris{PRN,4} -> Galileo (E)
%   ephemeris{PRN,5} -> QZSS (J)

% --- Initialize storage ---
ephemeris = cell(63, 5);

fid = fopen(filepath, 'r');
if fid == -1
    error('Cannot open file: %s', filepath);
end

cleanupObj = onCleanup(@() fclose(fid));

% --- Read header ---
while true
    line = fgetl(fid);
    if ~ischar(line)
        error('RINEX NAV header incomplete or END OF HEADER not found.');
    end
    if contains(line, 'END OF HEADER')
        break;
    end
end

% --- Read records ---
while ~feof(fid)
    line_raw = fgetl(fid);
    if ~ischar(line_raw)
        break;
    end
    if isempty(strtrim(line_raw))
        continue;
    end

    line = normalize_nav_line(line_raw);
    sys_char = upper(line(1));

    if ~ismember(sys_char, ['G', 'C', 'R', 'E', 'J'])
        continue;
    end

    try
        eph = struct();
        sys_idx = [];

        eph.PRN = str2double(line(2:3));
        if isnan(eph.PRN) || eph.PRN < 1 || eph.PRN > 63
            continue;
        end

        eph.Toc.Year = str2double(line(5:8));
        eph.Toc.Month = str2double(line(10:11));
        eph.Toc.Day = str2double(line(13:14));
        eph.Toc.Hour = str2double(line(16:17));
        eph.Toc.Minute = str2double(line(19:20));
        eph.Toc.Second = str2double(line(22:23));

        switch sys_char
            case 'G'
                eph.System = 'GPS';
                sys_idx = 1;

                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));

                for i = 1:7
                    line = read_nav_line(fid);
                    switch i
                        case 1
                            eph.IODE = str2double(line(5:23));
                            eph.Crs = str2double(line(24:42));
                            eph.Delta_n = str2double(line(43:61));
                            eph.M0 = str2double(line(62:80));
                        case 2
                            eph.Cuc = str2double(line(5:23));
                            eph.e = str2double(line(24:42));
                            eph.Cus = str2double(line(43:61));
                            eph.sqrtA = str2double(line(62:80));
                        case 3
                            eph.Toe = str2double(line(5:23));
                            eph.Cic = str2double(line(24:42));
                            eph.OMEGA0 = str2double(line(43:61));
                            eph.Cis = str2double(line(62:80));
                        case 4
                            eph.i0 = str2double(line(5:23));
                            eph.Crc = str2double(line(24:42));
                            eph.omega = str2double(line(43:61));
                            eph.OMEGA_DOT = str2double(line(62:80));
                        case 5
                            eph.IDOT = str2double(line(5:23));
                            eph.Codes_on_L2 = str2double(line(24:42));
                            eph.GPS_Week = str2double(line(43:61));
                            eph.L2_P_data_flag = str2double(line(62:80));
                        case 6
                            eph.SV_accuracy = str2double(line(5:23));
                            eph.SV_health = str2double(line(24:42));
                            eph.TGD = str2double(line(43:61));
                            eph.IODC = str2double(line(62:80));
                        case 7
                            eph.Transmission_time = str2double(line(5:23));
                            eph.Fit_interval = str2double(line(24:42));
                    end
                end

            case 'C'
                eph.System = 'BeiDou';
                sys_idx = 2;

                eph.A0 = str2double(line(24:42));
                eph.A1 = str2double(line(43:61));
                eph.A2 = str2double(line(62:80));

                for i = 1:7
                    line = read_nav_line(fid);
                    switch i
                        case 1
                            eph.IODE = str2double(line(5:23));
                            eph.Crs = str2double(line(24:42));
                            eph.Delta_n = str2double(line(43:61));
                            eph.M0 = str2double(line(62:80));
                        case 2
                            eph.Cuc = str2double(line(5:23));
                            eph.e = str2double(line(24:42));
                            eph.Cus = str2double(line(43:61));
                            eph.sqrtA = str2double(line(62:80));
                        case 3
                            eph.Toe = str2double(line(5:23));
                            eph.Cic = str2double(line(24:42));
                            eph.OMEGA0 = str2double(line(43:61));
                            eph.Cis = str2double(line(62:80));
                        case 4
                            eph.i0 = str2double(line(5:23));
                            eph.Crc = str2double(line(24:42));
                            eph.omega = str2double(line(43:61));
                            eph.OMEGA_DOT = str2double(line(62:80));
                        case 5
                            eph.IDOT = str2double(line(5:23));
                            eph.Spare1 = str2double(line(24:42));
                            eph.BDS_Week = str2double(line(43:61));
                            eph.Spare2 = str2double(line(62:80));
                        case 6
                            eph.SV_accuracy = str2double(line(5:23));
                            eph.SatH1 = str2double(line(24:42));
                            eph.TGD1 = str2double(line(43:61));
                            eph.TGD2 = str2double(line(62:80));
                        case 7
                            eph.Transmission_time = str2double(line(5:23));
                            eph.AODC = str2double(line(24:42));
                    end
                end

            case 'R'
                eph.System = 'GLONASS';
                sys_idx = 3;

                eph.tau_n = str2double(line(24:42));
                eph.gamma_n = str2double(line(43:61));
                eph.tk = str2double(line(62:80));

                for i = 1:3
                    line = read_nav_line(fid);
                    switch i
                        case 1
                            eph.pos_x = str2double(line(5:23));
                            eph.vel_x = str2double(line(24:42));
                            eph.acc_x = str2double(line(43:61));
                            eph.health = str2double(line(62:80));
                        case 2
                            eph.pos_y = str2double(line(5:23));
                            eph.vel_y = str2double(line(24:42));
                            eph.acc_y = str2double(line(43:61));
                            eph.freq_num = str2double(line(62:80));
                        case 3
                            eph.pos_z = str2double(line(5:23));
                            eph.vel_z = str2double(line(24:42));
                            eph.acc_z = str2double(line(43:61));
                            eph.age = str2double(line(62:80));
                    end
                end

            case 'E'
                eph.System = 'Galileo';
                sys_idx = 4;

                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));

                for i = 1:7
                    line = read_nav_line(fid);
                    switch i
                        case 1
                            eph.IODE = str2double(line(5:23));
                            eph.Crs = str2double(line(24:42));
                            eph.Delta_n = str2double(line(43:61));
                            eph.M0 = str2double(line(62:80));
                        case 2
                            eph.Cuc = str2double(line(5:23));
                            eph.e = str2double(line(24:42));
                            eph.Cus = str2double(line(43:61));
                            eph.sqrtA = str2double(line(62:80));
                        case 3
                            eph.Toe = str2double(line(5:23));
                            eph.Cic = str2double(line(24:42));
                            eph.OMEGA0 = str2double(line(43:61));
                            eph.Cis = str2double(line(62:80));
                        case 4
                            eph.i0 = str2double(line(5:23));
                            eph.Crc = str2double(line(24:42));
                            eph.omega = str2double(line(43:61));
                            eph.OMEGA_DOT = str2double(line(62:80));
                        case 5
                            eph.IDOT = str2double(line(5:23));
                            eph.Data_sources = str2double(line(24:42));
                            eph.GAL_Week = str2double(line(43:61));
                        case 6
                            eph.SV_accuracy = str2double(line(5:23));
                            eph.SV_health = str2double(line(24:42));
                            eph.TGD_E1E5a = str2double(line(43:61));
                            eph.TGD_E1E5b = str2double(line(62:80));
                        case 7
                            eph.Transmission_time = str2double(line(5:23));
                            eph.Spare = str2double(line(24:42));
                    end
                end

            case 'J'
                eph.System = 'QZSS';
                sys_idx = 5;

                eph.af0 = str2double(line(24:42));
                eph.af1 = str2double(line(43:61));
                eph.af2 = str2double(line(62:80));

                for i = 1:7
                    line = read_nav_line(fid);
                    switch i
                        case 1
                            eph.IODE = str2double(line(5:23));
                            eph.Crs = str2double(line(24:42));
                            eph.Delta_n = str2double(line(43:61));
                            eph.M0 = str2double(line(62:80));
                        case 2
                            eph.Cuc = str2double(line(5:23));
                            eph.e = str2double(line(24:42));
                            eph.Cus = str2double(line(43:61));
                            eph.sqrtA = str2double(line(62:80));
                        case 3
                            eph.Toe = str2double(line(5:23));
                            eph.Cic = str2double(line(24:42));
                            eph.OMEGA0 = str2double(line(43:61));
                            eph.Cis = str2double(line(62:80));
                        case 4
                            eph.i0 = str2double(line(5:23));
                            eph.Crc = str2double(line(24:42));
                            eph.omega = str2double(line(43:61));
                            eph.OMEGA_DOT = str2double(line(62:80));
                        case 5
                            eph.IDOT = str2double(line(5:23));
                            eph.Codes_on_L2 = str2double(line(24:42));
                            eph.GPS_Week = str2double(line(43:61));
                        case 6
                            eph.SV_accuracy = str2double(line(5:23));
                            eph.SV_health = str2double(line(24:42));
                            eph.TGD = str2double(line(43:61));
                            eph.IODC = str2double(line(62:80));
                        case 7
                            eph.Transmission_time = str2double(line(5:23));
                            eph.Fit_interval = str2double(line(24:42));
                    end
                end
        end

        if ~isempty(sys_idx)
            if isempty(ephemeris{eph.PRN, sys_idx})
                ephemeris{eph.PRN, sys_idx} = eph;
            else
                ephemeris{eph.PRN, sys_idx}(end+1) = eph;
            end
        end

    catch ME
        sat_id = strtrim(line(1:3));
        if isempty(sat_id)
            sat_id = 'N/A';
        end
        fprintf('Warning: failed to parse satellite %s (system %s): %s. Skip this record.\n', ...
            sat_id, sys_char, ME.message);

        if strcmp(ME.identifier, 'PARSE_RINEX_NAV_MULTI_GNSS:UnexpectedEOF') || feof(fid)
            break;
        end
    end
end

fprintf('Done parsing multi-GNSS nav file (G=1, C=2, R=3, E=4, J=5).\n');
end

function line = read_nav_line(fid)
line = fgetl(fid);
if ~ischar(line)
    error('PARSE_RINEX_NAV_MULTI_GNSS:UnexpectedEOF', ...
        'Unexpected EOF while reading NAV record.');
end
line = normalize_nav_line(line);
end

function line = normalize_nav_line(line)
line = strrep(line, 'D', 'E');
if length(line) < 80
    line = [line, repmat(' ', 1, 80 - length(line))];
end
end
