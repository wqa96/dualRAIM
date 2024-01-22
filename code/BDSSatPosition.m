function [satdata] = BDSSatPosition(navdata, t, xyz_type)
% ������������navdata����tʱ�̣�UTC��������xyz_type����ϵ������ϵI����ĵع�ϵECEF���µ�λ��
% navdata = struct('prn', [], 'time', [], 'a012', [], ...
%         'aode', [], 'crs', [], 'dn', [], 'm0', [], ...
%         'cuc', [], 'e', [], 'cus', [], 'sqrtA', [], ...
%         'toe', [], 'cic', [], 'Omega0', [], 'cis', [], ...
%         'i0', [], 'crc', [], 'omega', [], 'dOmega', [], ...
%         'di', [], 'sp1', [], 'week', [], 'sp2', [], ...
%         'svaccuracy', [], 'sathl', [], 'tgd1', [], 'tgd2', [], ...
%         'transmission', [], 'aodc', [], 'sp3', [], 'sp4', []);
    satdata = struct('prn', [], ...
        'x', [], 'y', [], 'z', [], 'svaccuracy' , [], ...
        'vx', [], 'vy', [], 'vz', [], ...
        't', [], 'toe', [], 'week', [], ...
        'sathl', [], 'transmission', [], ...
        'aode', [], 'aodc', [], ...
        'tgd1', [], 'tgd2', [], 'a012', []);
    BDSLeapSeconds = 4;
    mu = 3.986005e14;
    omegae = 7.2921151467e-5;
    t1 = datetime([2006 1 1 0 0 0]);
    t2BDT = mod(seconds(datetime(t)-t1) + BDSLeapSeconds, 604800);
    tk = t2BDT - navdata.toe;
    ind = find((tk >= 0) & (tk < 3600));
    %����ɸѡ
    prn = str2num(navdata.prn(ind, 2:3));
    tk = tk(ind);
    toe = navdata.toe(ind);
    sqrtA = navdata.sqrtA(ind);
    A = sqrtA.^2;
    dn = navdata.dn(ind);
    m0 = navdata.m0(ind);
    e = navdata.e(ind);
    omega = navdata.omega(ind);
    cus = navdata.cus(ind);
    cuc = navdata.cuc(ind);
    crs = navdata.crs(ind);
    crc = navdata.crc(ind);
    cis = navdata.cis(ind);
    cic = navdata.cic(ind);
    i0 = navdata.i0(ind);
    di = navdata.di(ind);
    Omega0 = navdata.Omega0(ind);
    dOmega = navdata.dOmega(ind);
    %ƽ�����ٶ�
    n0 = sqrt(mu./sqrtA.^6);
    n = n0 + dn;
    %ƽ�����
    mk = m0 + n.*tk;
    dmk = n;
    %ƫ�����
    ek0 = 0;
    ek = mk;
    cnt = 0;
    while(all(abs(ek-ek0) > 1e-15))
        ek0 = ek;
        ek = mk + e.*sin(ek);
        cnt = cnt + 1;
        if(cnt > 1e3)
            error('ƫ����Ǽ����޷�����');
        end
    end
    dek = dmk ./ (1-e.*cos(ek));
    %������
    tmp0 = sqrt(1 - e.^2);
    tmp1 = 1 - e.*cos(ek);
    v1 = (cos(ek) - e)./tmp1;
    v2 = (tmp0.*sin(ek))./tmp1;
    vk = atan2(v2, v1);
    dvk = tmp0.*dek./tmp1;
    %�������
    phik = vk + omega;
    dphik = dvk;
    %�㶯У����
    sinPhik = sin(2*phik);
    cosPhik = cos(2*phik);
    duk = cus.*sinPhik + cuc.*cosPhik;
    drk = crs.*sinPhik + crc.*cosPhik;
    dik = cis.*sinPhik + cic.*cosPhik;
    dduk = 2*dphik.*(cus.*cosPhik - cuc.*sinPhik);
    ddrk = 2*dphik.*(crs.*cosPhik - crc.*sinPhik);
    ddik = 2*dphik.*(cis.*cosPhik - cic.*sinPhik);
    %������Ǿ�
    uk = phik + duk;
    duk0 = dphik + dduk;
    %����ʸ������
    rk = A.*tmp1 + drk;
    drk0 = A.*e.*dek.*sin(ek) + ddrk;
    %������
    ik = i0 + di.*tk + dik;
    dik0 = di + ddik;
    %�����ڹ��ƽ�������
    xk0 = rk.*cos(uk);
    yk0 = rk.*sin(uk);
    vxk0 = drk0.*cos(uk) - rk.*duk0.*sin(uk);
    vyk0 = drk0.*sin(uk) + rk.*duk0.*cos(uk);

    xk = zeros(size(tk));
    yk = zeros(size(tk));
    zk = zeros(size(tk));
    vxk = zeros(size(tk));
    vyk = zeros(size(tk));
    vzk = zeros(size(tk));
    if strcmp(xyz_type, 'ECEF')
        %% MEO/IGSO��������
        ind0 = (prn <= 5) | (prn >= 59);%GEO
        %������ྭ
        Omegak = Omega0 + (dOmega-omegae).*tk - omegae*toe;
        dOmegak = dOmega - omegae;
        %GEO������ྭ
        Omegak(ind0) = Omegak(ind0) + omegae.*tk(ind0);
        dOmegak(ind0) = dOmegak(ind0) + omegae;
        %ECEF
        xk = xk0.*cos(Omegak) - yk0.*sin(Omegak).*cos(ik);
        yk = xk0.*sin(Omegak) + yk0.*cos(Omegak).*cos(ik);
        zk = yk0.*sin(ik);
        vxk = -yk.*dOmegak + vxk0.*cos(Omegak) - ...
            (vyk0.*cos(ik) - zk.*dik0).*sin(Omegak);
        vyk = xk.*dOmegak + vxk0.*sin(Omegak) + ...
            (vyk0.*cos(ik) - zk.*dik0).*cos(Omegak);
        vzk = vyk0.*sin(ik) + yk0.*dik0.*cos(ik);
        %% GEO��������
        ind1 = find(ind0);
        Rx = [1, 0, 0;
            0, cos(-5/180*pi), sin(-5/180*pi);
            0, -sin(-5/180*pi), cos(-5/180*pi)];
        xyzk = zeros(3, length(ind1));
        vxyzk = zeros(3, length(ind1));
        for n = 1:length(ind1)
            Rz = [cos(omegae*tk(ind1(n))), sin(omegae*tk(ind1(n))), 0;
               -sin(omegae*tk(ind1(n))), cos(omegae*tk(ind1(n))), 0;
               0, 0, 1];
            xyzk(:, n) = Rz*Rx*[xk(n); yk(n); zk(n)];
            vxyzk(:, n) = Rz*Rx*[vxk(n); vyk(n); vzk(n)] + ...
                [xyzk(2, n)*omegae; -xyzk(1, n)*omegae; 0];
        end
        %ECEF
        xk(ind0) = xyzk(1, :);
        yk(ind0) = xyzk(2, :);
        zk(ind0) = xyzk(3, :);
        vxk(ind0) = vxyzk(1, :);
        vyk(ind0) = vxyzk(2, :);
        vzk(ind0) = vxyzk(3, :);
    elseif strcmp(xyz_type, 'I')
        %% MEO/IGSO��������
        %������ྭ
        Omegak = Omega0 + dOmega.*tk;
        %����ϵ
        xk = xk0.*cos(Omegak) - yk0.*sin(Omegak).*cos(ik);
        yk = xk0.*sin(Omegak) + yk0.*cos(Omegak).*cos(ik);
        zk = yk0.*sin(ik);
        %% GEO��������
        ind0 = (prn <= 5) | (prn >= 59);%GEO
        Rx = [1, 0, 0;
            0, cos(-5/180*pi), sin(-5/180*pi);
            0, -sin(-5/180*pi), cos(-5/180*pi)];
%         Rx = [1, 0, 0;
%             0, 1, 0;
%             0, 0, 1];
        xyzk = Rx*[xk(ind0)'; yk(ind0)'; zk(ind0)'];
        %����ϵ
        xk(ind0) = xyzk(1, :);
        yk(ind0) = xyzk(2, :);
        zk(ind0) = xyzk(3, :);
    else
        disp('��������ϵ����xyz_typeӦ��Ϊ�����ĵع�ϵECEF�����ϵI');
        return;
    end
    satdata.prn = navdata.prn(ind, :);
    satdata.x = xk;
    satdata.y = yk;
    satdata.z = zk;
    satdata.t = t2BDT;
    satdata.vx = vxk;
    satdata.vy = vyk;
    satdata.vz = vzk;
    satdata.svaccuracy = navdata.svaccuracy(ind);
    satdata.toe = navdata.toe(ind);
    satdata.week = navdata.week(ind);
    satdata.sathl = navdata.sathl(ind);
    satdata.transmission = navdata.transmission(ind);
    satdata.aode = navdata.aode(ind);
    satdata.aodc = navdata.aodc(ind);
    satdata.tgd1 = navdata.tgd1(ind);
    satdata.tgd2 = navdata.tgd2(ind);
    satdata.a012 = navdata.a012(ind);
end