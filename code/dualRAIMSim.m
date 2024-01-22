% �޳�������2���쳣α����㷨����
clear;
% close all;
tic;
addpath ..\data;
c = 299792458;              % ���٣�m/s��
%% ����ģʽ
N = 10000;                  %�������simulation counts
satnum = 12; %satellite number9, 12, 15, 20
alarmRateInd = 3;
spSatNmax = 2;
prangeBias = [100, 1e4];
%% ��ȡ�ɼ���������get satellite informations
navdata = BDSeph('hour0010.20b', '3.03');   %��ȡĳһ��ı��������ļ�rinex3.03 read an rinex3.03 file
PRN = unique(navdata.prn, 'rows');          %��ȡ���е����Ǳ��CXX get the prn number
PRNNUM = str2num(PRN(:, 2:3));              %PRN���תΪ����XX

refpLLA = [108.930546733333, 34.198421205, 1413.8541];   %�ο��㾭γ��the reference LLA position(long, lat, alt)
[refpX, refpY, refpZ] = LLA2ECEF0(refpLLA(1), refpLLA(2), refpLLA(3));   %�ο���ECEF the reference ECEF position(X, Y, Z)

t = [2020 1 1 1 1 0];                       %����ʱ��reference time

satdata = BDSSatPosition(navdata, t, 'ECEF');   %����ʱ�̵��������� the satellites' data at reference time
sat2ref = [satdata.x - refpX, ...
    satdata.y - refpY, ...
    satdata.z - refpZ];                     %����λ�ù��ڲο�λ�õ�ʸ��
r = sqrt(sum(sat2ref.^2, 2));               %ʸ������vector length
satind = ((sat2ref * [refpX; refpY; refpZ]) > 0) & (satdata.sathl == 0);%�ο��㴦�ɼ�����satellites in view
% MEO_IGSO = (PRNNUM > 5) & (PRNNUM < 59);    %MEO/IGSO����
if satnum == 9
    satind([4, 5, 9, 14, 15, 19, 23, 30, 32, 34, 35]) = false;
elseif satnum == 12
    satind([14, 15, 19, 23, 30, 32, 34, 35]) = false;
elseif satnum == 15
    satind([23, 30, 32, 34, 35]) = false;
end
prnnum = PRNNUM(satind);                    %�ο��㴦�ɼ����Ǳ��satellites' prn in view
satlabels = mat2cell(satdata.prn(satind, :), ones(sum(satind), 1), 3);
%% �������/���ͳ����/����
%�������variance of pseudorange measurement noise
sigmaURE = 4;
%���ͳ������ֵthe SSE threshold of diffetent alarm rate
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
pos_r = [refpX; refpY; refpZ];%real position
%��������
satdata = BDSSatPosition(navdata, [2020 1 1 1 1 0], 'ECEF');    %����ʱ�̵���������
satpos = [satdata.x(satind), satdata.y(satind), satdata.z(satind)]';
%��ʵ����
r = mynorm(satpos - pos_r, 2, 1)';%satellites' geometric distance
%������ת����add earth rotation error
earth_rot = 7.292115e-5 * (satpos(1, :) * pos_r(2) - satpos(2, :) * pos_r(1))' / c;
r = r + earth_rot;
%% MontoCarlo����
%��ƭ�������all combinations of satellites (the satellite number no less than 5)
spSatInd = logical(dec2bin(1:2^satnum-2, satnum) - '0');
spSatNum = sum(spSatInd, 2);
[spSatNum, ind] = sort(spSatNum, 'ascend');
spSatInd = spSatInd(ind, :);
spSatL2 = spSatNum <= spSatNmax;
spSatInd = spSatInd(spSatL2, :);
spSatNum = spSatNum(spSatL2);
% ������progress bar
tElapsed = 0;
spSatIdNum = length(spSatNum);
simTotNum = N*spSatIdNum;
simCurNum = 0;
simRemNum = simTotNum;
waitbarWin = waitbar(0, '��ʼ��');
%�����ж�the result of the algorithm
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
        % ������progress bar
        simCurNum = N*(spSatId-1)+n;
        simRemNum = simTotNum - simCurNum;
        waitbar(simCurNum/simTotNum, waitbarWin, ...
            ['�������: ', num2str(simCurNum/simTotNum*100, '%.2f'), ' % (', ...
            num2str(simCurNum), '/', num2str(simTotNum), ')', ...
            'ʣ��ʱ��:', num2str(remHour), ':', num2str(remMin, '%02d'), ':', num2str(remSec, '%02d')]);
        %% α��۲�������
        %���ջ�ʱ�����time bias
        dtr = randn()*1e-4;
        %����α��pseudorange
        pr = r;
        pr(spInd) = pr(spInd) + prangeBias(1) + (prangeBias(2)-prangeBias(1))*rand(spn, 1);
        URE = sigmaURE*randn(satnum, 1);
        pr = pr + c*dtr + URE;
        %��λ����pvt
        pos = [pos_r; 0]';
%         [pos, rtmp, pErr, G, H, cnt] = pvtpos0(pr(~spInd), satpos(:, ~spInd)', pos, 1e-10, 10);
        [pos, rtmp, pErr, G, H, cnt] = pvtpos0(pr, satpos', pos, 1e-10, 10);
        %% dual-RAIM
        spDet = false(satnum, 1);   %��ƭ�����fault satellite flag
        spDet0 = false(satnum, 1);   %��ƭ�����fault satellite flag
        if pErr'*pErr > SSE(alarmRateInd, satnum-4)
            % ��⵽��ƭdetect the existance of failure
            S = eye(satnum) - G*H*G';
            Fi = pErr.^2./diag(S);
            [Fimax, FiInd] = max(Fi);
            SSEi = pErr'*pErr - Fimax;
            if (SSEi < SSE(alarmRateInd, satnum-5))
                % ��⵽�����ڵ����쳣α��find only one fault
                spDet(FiInd) = true;
            else
                % ��⵽����쳣α��find more than one fault
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
                    % ��⵽2���쳣α��find two faults
                    spDet([FijInd1, FijInd2]) = true;
                else
                    % ��⵽����2���쳣α��find more than two faults
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
            % δ��⵽��ƭno faults
            detect(spSatId, n) = 0;
        elseif all(spDet == spInd)
            % �ɹ������������ƭfind all faults
            detect(spSatId, n) = 1;
        elseif all(spDet)
            % ��Ϊ���ڶ���2���쳣����think there are more than two faults
            detect(spSatId, n) = 2;
        elseif any(spDet(~spInd)) && any(spInd(~spDet))
            % ��������������Ϊ�쳣���ǲ������쳣����û�м�����
            % a right measurement detect as a fault and there is a fault
            % not detected
            detect(spSatId, n) = 3;
        elseif any(spDet(~spInd))
            % ��������������Ϊ�쳣����
            % a right measurement detect as a fault
            detect(spSatId, n) = 4;
        elseif any(spInd(~spDet))
            % ���쳣����û�м�����
            % there is a fault not detected
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
%% ���ݱ���save data
t = datetime('now');
tMDHM = reshape(num2str([t.Month, t.Day, t.Hour, t.Minute]', '%02d')', 1, []);
fullpath = mfilename('fullpath');
[~, name] = fileparts(fullpath);
savefilename = ['..\data\', name, ...
    '_', num2str([satnum, spSatNmax], '%02d'), ...
    '_', num2str(alarmRateInd), ...
    '_', [num2str(prangeBias(2)), 'm'], ...
    '_', tMDHM, '.mat'];
save(savefilename);
%% ����SSE����
figInd = 1;

failDet = zeros(spSatNmax, 1);
falseDet = zeros(spSatNmax, 4);
falseNums = zeros(3, 4, spSatNmax);
trueDet = zeros(spSatNmax, 1);
falseProb = 0;
failProb = 0;
for n = 1:spSatNmax
    spnind = (spSatNum == n);
    dettmp = detect(spnind, :);
    trueDet(n) = sum(dettmp(:) == 1) / (N*sum(spnind));
    failDet(n) = sum(dettmp(:) == 0) / (N*sum(spnind));
    failProb = failProb + sum(dettmp(:) == 0);
    falseDet(n, :) = sum(dettmp(:) == (2:5), 1) / (N*sum(spnind));
    for nn = 1:3
        indtmp = dettmp == nn+2;
        datatmp = falseNum(spnind, :);
        falseNums(nn, :, n) = sum(datatmp(indtmp) == (1:4), 1) / (N*sum(spnind));
    end
    falseProb = falseProb + sum(dettmp(:) > 1);
end
falseProb = falseProb / (N*spSatIdNum);
failProb = failProb / (N*spSatIdNum);
fprintf('ƽ���ɹ����ʣ�%.7f\tƽ�������ʣ�%.7f\tƽ��ʧ�ܸ��ʣ�%.7f\n', 1-falseProb-failProb ,falseProb, failProb);

figure(figInd); figInd = figInd + 1;
bar(1:spSatNmax, [trueDet, sum(falseDet, 2), failDet], 'stacked'); 
grid on; xlim([0, spSatNmax+1]); ylim([0, 1]);
legend('success', 'false', 'fail', 'Location', 'bestoutside');
xlabel('number of spoofing signals'), ylabel('probability');


figure(figInd); figInd = figInd + 1;
bar(1:spSatNmax, falseDet, 'stacked'); grid on;
xlim([0, spSatNmax+1]); title('false detection');
legend('SatM2', 'False&Miss', 'False', 'Miss', ...
    'Location','bestoutside');
xlabel('number of spoofing signals'), ylabel('probability');
%% ��ȷ����SSE����
failDet = zeros(spSatNmax, 1);
falseDet = zeros(spSatNmax, 4);
trueDet = zeros(spSatNmax, 1);
falseProb = 0;
failProb = 0;
for n = 1:spSatNmax
    spnind = (spSatNum == n);
    dettmp = detect0(spnind, :);
    trueDet(n) = sum(dettmp(:) == 1) / (N*sum(spnind));
    failDet(n) = sum(dettmp(:) == 0) / (N*sum(spnind));
    failProb = failProb + sum(dettmp(:) == 0);
    falseDet(n, :) = sum(dettmp(:) == (2:5), 1) / (N*sum(spnind));
    falseProb = falseProb + sum(dettmp(:) > 1);
end
falseProb = falseProb / (N*spSatIdNum);
failProb = failProb / (N*spSatIdNum);
fprintf('ƽ�������ʣ�%.4f\tƽ��ʧ�ܸ��ʣ�%.4f\n', falseProb, failProb);

figure(figInd); figInd = figInd + 1;
bar(1:spSatNmax, [trueDet, sum(falseDet, 2), failDet], 'stacked'); 
grid on; xlim([0, spSatNmax+1]); ylim([0, 1]);
legend('success', 'false', 'fail', 'Location', 'bestoutside');
xlabel('number of spoofing signals'), ylabel('probability');


figure(figInd); figInd = figInd + 1;
bar(1:spSatNmax, falseDet, 'stacked'); grid on;
xlim([0, spSatNmax+1]); title('false detection');
legend('SatM2', 'False&Miss', 'False', 'Miss', ...
    'Location','bestoutside');
xlabel('number of spoofing signals'), ylabel('probability');
