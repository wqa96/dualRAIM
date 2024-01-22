function [navdata] = BDSeph(filename, rinex_type)
% ½âÎö±±¶·ÐÇÀúrinex3.03
    navdata = struct('prn', [], 'time', [], 'a012', [], ...
            'aode', [], 'crs', [], 'dn', [], 'm0', [], ...
            'cuc', [], 'e', [], 'cus', [], 'sqrtA', [], ...
            'toe', [], 'cic', [], 'Omega0', [], 'cis', [], ...
            'i0', [], 'crc', [], 'omega', [], 'dOmega', [], ...
            'di', [], 'sp1', [], 'week', [], 'sp2', [], ...
            'svaccuracy' ,[], 'sathl', [], 'tgd1', [], 'tgd2', [], ...
            'transmission', [], 'aodc', [], 'sp3', [], 'sp4', []);
    if ~exist(filename, 'file')
        disp([filename, ' not found']);
        return;
    end
    if strcmp(rinex_type, '3.03')
        fileID = fopen(filename, 'r');
        while(1)
            tline = fgetl(fileID);
            if(tline == -1)
                disp('END OF HEADER not found');
                fclose(fileID);
                return;
            end
            if contains(tline, 'END OF HEADER')
                break;
            end
        end
        C = textscan(fileID, ['%3c %d %d %d %d %d %d %f %f %f', ...
            ' %f %f %f %f', ' %f %f %f %f', ' %f %f %f %f', ...
            ' %f %f %f %f', ' %f %f %f %f', ' %f %f %f %f', ...
            ' %f %f'], 'Delimiter', '\n');
        navdata.prn = C{1};
        navdata.time = [C{2}, C{3}, C{4}, C{5}, C{6}, C{7}];
        navdata.a012 = [C{8}, C{9}, C{10}];
        navdata.aode = C{11};
        navdata.crs = C{12};
        navdata.dn = C{13};
        navdata.m0 = C{14};
        navdata.cuc = C{15};
        navdata.e = C{16};
        navdata.cus = C{17};
        navdata.sqrtA = C{18};
        navdata.toe = C{19};
        navdata.cic = C{20};
        navdata.Omega0 = C{21};
        navdata.cis = C{22};
        navdata.i0 = C{23};
        navdata.crc = C{24};
        navdata.omega = C{25};
        navdata.dOmega = C{26};
        navdata.di = C{27};
        navdata.sp1 = C{28};
        navdata.week = C{29};
        navdata.sp2 = C{30};
        navdata.svaccuracy = C{31};
        navdata.sathl = C{32};
        navdata.tgd1 = C{33};
        navdata.tgd2 = C{34};
        navdata.transmission = C{35};
        navdata.aodc = C{36};
        navdata.sp3 = nan(length(C{1}), 1);
        navdata.sp4 = nan(length(C{1}), 1);
        fclose(fileID);
    else
        disp(['rinex type:', rinex_type, ' not support']);
    end
end