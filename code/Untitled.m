clear;
addpath data;

N = 10000;
satnum = 12;
spSatNmax = 2;
prangeBias = [100, 1e4;
    50, 100;
    30, 50;
    10, 30];
% parfor alarmRateInd = 1:6
%     dualRAIMSimulation(N, satnum, alarmRateInd, spSatNmax, prangeBias);
% end
parfor prangeBiasInd = 1:4
    dualRAIMSimulation(N, satnum, 2, spSatNmax, prangeBias(prangeBiasInd, :));
end
%%
function dualRAIMSimulation(N, satnum, alarmRateInd, spSatNmax, prangeBias)
    c = 299792458;              % ���٣�m/s��
    %% ��ȡ�ɼ���������
    navdata = BDSeph('hour0010.20b', '3.03');   %��ȡĳһ��ı��������ļ�rinex3.03
    PRN = unique(navdata.prn, 'rows');          %��ȡ���е����Ǳ��CXX
    PRNNUM = str2num(PRN(:, 2:3));              %PRN���תΪ����XX

    refpLLA = [108.930546733333, 34.198421205, 1413.8541];   %�ο��㾭γ��
    [refpX, refpY, refpZ] = LLA2ECEF0(refpLLA(1), refpLLA(2), refpLLA(3));   %�ο���ECEF

    t = [2020 1 1 1 1 0];                       %����ʱ��

    satdata = BDSSatPosition(navdata, t, 'ECEF');   %����ʱ�̵���������
    sat2ref = [satdata.x - refpX, ...
        satdata.y - refpY, ...
        satdata.z - refpZ];                     %����λ�ù��ڲο�λ�õ�ʸ��
    r = sqrt(sum(sat2ref.^2, 2));               %ʸ������
    satind = ((sat2ref * [refpX; refpY; refpZ]) > 0) & (satdata.sathl == 0);%�ο��㴦�ɼ�����
    % MEO_IGSO = (PRNNUM > 5) & (PRNNUM < 59);    %MEO/IGSO����
    if satnum == 9
        satind([4, 5, 9, 14, 15, 19, 23, 30, 32, 34, 35]) = false;
    elseif satnum == 12
        satind([14, 15, 19, 23, 30, 32, 34, 35]) = false;
    elseif satnum == 15
        satind([23, 30, 32, 34, 35]) = false;
    end
    prnnum = PRNNUM(satind);                    %�ο��㴦�ɼ����Ǳ��
    satlabels = mat2cell(satdata.prn(satind, :), ones(sum(satind), 1), 3);
    %% �������/���ͳ����/����
    %�������
    sigmaURE = 4;
    %���ͳ������ֵ
    alarmRate = [0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001];
    SSE = zeros(length(alarmRate), 17);
    alpha = 0:0.01:60;
    for m = 1:length(alarmRate)
        for n = 1:17
            [B, I] = min(abs(cdf('Chisquare', alpha, n) - (1-alarmRate(m))));
            SSE(m, n) = alpha(I)*sigmaURE^2;
        end
    end
    %���ͳ���������ж���1��������
    SPAlarm = @(x, n, m) x'*x > SSE(m, n);
    %% �켣/����/��ʵ����
    %λ��
    pos_r = [refpX; refpY; refpZ];
    %��������
    satdata = BDSSatPosition(navdata, [2020 1 1 1 1 0], 'ECEF');    %����ʱ�̵���������
    satpos = [satdata.x(satind), satdata.y(satind), satdata.z(satind)]';
    %��ʵ����
    r = mynorm(satpos - pos_r, 2, 1)';
    %������ת����
    earth_rot = 7.292115e-5 * (satpos(1, :) * pos_r(2) - satpos(2, :) * pos_r(1))' / c;
    %��ʵα��
    r = r + earth_rot;
    %% MontoCarlo����
    %��ƭ�������
    spSatInd = logical(dec2bin(1:2^satnum-2, satnum) - '0');
    spSatNum = sum(spSatInd, 2);
    [spSatNum, ind] = sort(spSatNum, 'ascend');
    spSatInd = spSatInd(ind, :);
    spSatL2 = spSatNum <= spSatNmax;
    spSatInd = spSatInd(spSatL2, :);
    spSatNum = spSatNum(spSatL2);
    % ������
    tElapsed = 0;
    spSatIdNum = length(spSatNum);
    simTotNum = N*spSatIdNum;
    simCurNum = 0;
    simRemNum = simTotNum;
    waitbarWin = waitbar(0, '��ʼ��');
    %�����ж�
    detect = zeros(spSatIdNum, N);       %�������Ƿ���ȷ����
    detect0 = zeros(spSatIdNum, N);       %�������Ƿ���ȷ����
    falseNum = zeros(spSatIdNum, N);
    falseNum0 = zeros(spSatIdNum, N);
    tic;
    for spSatId = 1:spSatIdNum
        n = 1;
        spInd = spSatInd(spSatId, :)';
        spn = sum(spInd);

        simSpeed = simCurNum/tElapsed;
        remSec = simRemNum/simSpeed;
        remHour = floor(remSec/3600);
        remMin = floor(mod(remSec, 3600)/60);
        remSec = floor(mod(remSec, 60));
        while n <= N
            % ������
            simCurNum = N*(spSatId-1)+n;
            simRemNum = simTotNum - simCurNum;
            waitbar(simCurNum/simTotNum, waitbarWin, ...
                ['�������: ', num2str(simCurNum/simTotNum*100, '%.2f'), ' % (', ...
                num2str(simCurNum), '/', num2str(simTotNum), ')', ...
                'ʣ��ʱ��:', num2str(remHour), ':', num2str(remMin, '%02d'), ':', num2str(remSec, '%02d')]);
            %% α��۲�������
            %���ջ�ʱ�����
            dtr = randn()*1e-4;
            %����α��
            pr = r;
            pr(spInd) = pr(spInd) + prangeBias(1) + (prangeBias(2)-prangeBias(1))*rand(spn, 1);
            URE = sigmaURE*randn(satnum, 1);
            pr = pr + c*dtr + URE;
            %��λ����
            pos = [pos_r; 0]';
    %         [pos, rtmp, pErr, G, H, cnt] = pvtpos0(pr(~spInd), satpos(:, ~spInd)', pos, 1e-10, 10);
            [pos, rtmp, pErr, G, H, cnt] = pvtpos0(pr, satpos', pos, 1e-10, 10);
            %% dual-RAIM
            spDet = false(satnum, 1);   %��ƭ�����
            spDet0 = false(satnum, 1);   %��ƭ�����
            if pErr'*pErr > SSE(alarmRateInd, satnum-4)
                % ��⵽��ƭ
                S = eye(satnum) - G*H*G';
                Fi = pErr.^2./diag(S);
                [Fimax, FiInd] = max(Fi);
                SSEi = pErr'*pErr - Fimax;
                if (SSEi < SSE(alarmRateInd, satnum-5))
                    % ��⵽�����ڵ����쳣α��
                    spDet(FiInd) = true;
                else
                    % ��⵽����쳣α��
                    Fij = zeros(satnum, satnum);
                    for ki = 1:satnum
                        for kj = ki+1:satnum
                            Fij(ki, kj) = (S(ki,ki)*pErr(kj)^2 - 2*S(ki,kj)*pErr(ki)*pErr(kj) + S(kj,kj)*pErr(ki)^2)/(S(ki,ki)*S(kj,kj)-S(ki,kj)^2);
                        end
                    end
                    [Fijmax, FijInd] = max(Fij(:));
                    [FijInd1, FijInd2] = ind2sub(size(Fij),FijInd);
                    SSEij = pErr'*pErr - Fijmax;
                    if (SSEij < SSE(alarmRateInd, satnum-6))
                        % ��⵽2���쳣α��
                        spDet([FijInd1, FijInd2]) = true;
                    else
                        % ��⵽����2���쳣α��
                        spDet(:) = true;
                    end
                end
                %% ��ȷ�����޳����Ǻ��SSE
                [Fimax, FiInd] = max(Fi);
                parttmp = true(satnum, 1);
                parttmp(FiInd) = false;
                [pos_tmp, rtmp, pErrtmp, G, H] = pvtpos0(pr(parttmp), satpos(:, parttmp)', pos, 1e-10, 10);
                SSEi = pErrtmp'*pErrtmp;
                if (SSEi < SSE(alarmRateInd, satnum-5))
                    spDet0(FiInd) = true;
                else
                    Fij = zeros(satnum, satnum);
                    for ki = 1:satnum
                        for kj = ki+1:satnum
                            Fij(ki, kj) = (S(ki,ki)*pErr(kj)^2 - 2*S(ki,kj)*pErr(ki)*pErr(kj) + S(kj,kj)*pErr(ki)^2)/(S(ki,ki)*S(kj,kj)-S(ki,kj)^2);
                        end
                    end
                    [Fijmax, FijInd] = max(Fij(:));
                    [FijInd1, FijInd2] = ind2sub(size(Fij),FijInd);
                    parttmp = true(satnum, 1);
                    parttmp([FijInd1, FijInd2]) = false;
                    [pos_tmp, rtmp, pErrtmp, G, H] = pvtpos0(pr(parttmp), satpos(:, parttmp)', pos, 1e-10, 10);
                    SSEij = pErrtmp'*pErrtmp;
                    if (SSEij < SSE(alarmRateInd, satnum-6))
                        spDet0([FijInd1, FijInd2]) = true;
                    else
                        spDet0(:) = true;
                    end
                end
            end
            %% �����������
            falseNum(spSatId, n) = sum(xor(spDet, spInd));
            if all(~spDet)
                % δ��⵽��ƭ
                detect(spSatId, n) = 0;
            elseif all(spDet == spInd)
                % �ɹ������������ƭ
                detect(spSatId, n) = 1;
            elseif all(spDet)
                % ��Ϊ���ڶ���2���쳣����
                detect(spSatId, n) = 2;
            elseif any(spDet(~spInd)) && any(spInd(~spDet))
                % ��������������Ϊ�쳣���ǲ������쳣����û�м�����
                detect(spSatId, n) = 3;
            elseif any(spDet(~spInd))
                % ��������������Ϊ�쳣����
                detect(spSatId, n) = 4;
            elseif any(spInd(~spDet))
                % ���쳣����û�м�����
                detect(spSatId, n) = 5;
            end
            %% �����������1
            falseNum0(spSatId, n) = sum(xor(spDet0, spInd));
            if all(~spDet0)
                % δ��⵽��ƭ
                detect0(spSatId, n) = 0;
            elseif all(spDet0 == spInd)
                % �ɹ������������ƭ
                detect0(spSatId, n) = 1;
            elseif all(spDet0)
                % ��Ϊ���ڶ���2���쳣����
                detect0(spSatId, n) = 2;
            elseif any(spDet0(~spInd)) && any(spInd(~spDet0))
                % ��������������Ϊ�쳣���ǲ������쳣����û�м�����
                detect0(spSatId, n) = 3;
            elseif any(spDet0(~spInd))
                % ��������������Ϊ�쳣����
                detect0(spSatId, n) = 4;
            elseif any(spInd(~spDet0))
                % ���쳣����û�м�����
                detect0(spSatId, n) = 5;
            end
            n = n + 1;
        end
        tElapsed = toc;
    end
    close(waitbarWin);
    toc;
    %% ���ݱ���
    t = datetime('now');
    tMDHM = reshape(num2str([t.Month, t.Day, t.Hour, t.Minute]', '%02d')', 1, []);
    fullpath = mfilename('fullpath');
    [~, name] = fileparts(fullpath);
    savefilename = ['.\data\', name, ...
        '_', num2str([satnum, spSatNmax], '%02d'), ...
        '_', num2str(alarmRateInd), ...
        '_', [num2str(prangeBias(2)), 'm'], ...
        '_', tMDHM, '.mat'];
    save(savefilename);
end