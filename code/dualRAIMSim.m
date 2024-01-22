% 剔除不超过2个异常伪距的算法仿真
clear;
% close all;
tic;
addpath ..\data;
c = 299792458;              % 光速（m/s）
%% 仿真模式
N = 10000;                  %仿真次数simulation counts
satnum = 12; %satellite number9, 12, 15, 20
alarmRateInd = 3;
spSatNmax = 2;
prangeBias = [100, 1e4];
%% 获取可见卫星数据get satellite informations
navdata = BDSeph('hour0010.20b', '3.03');   %读取某一天的北斗星历文件rinex3.03 read an rinex3.03 file
PRN = unique(navdata.prn, 'rows');          %获取所有的卫星编号CXX get the prn number
PRNNUM = str2num(PRN(:, 2:3));              %PRN编号转为数字XX

refpLLA = [108.930546733333, 34.198421205, 1413.8541];   %参考点经纬高the reference LLA position(long, lat, alt)
[refpX, refpY, refpZ] = LLA2ECEF0(refpLLA(1), refpLLA(2), refpLLA(3));   %参考点ECEF the reference ECEF position(X, Y, Z)

t = [2020 1 1 1 1 0];                       %仿真时刻reference time

satdata = BDSSatPosition(navdata, t, 'ECEF');   %仿真时刻的卫星数据 the satellites' data at reference time
sat2ref = [satdata.x - refpX, ...
    satdata.y - refpY, ...
    satdata.z - refpZ];                     %卫星位置关于参考位置的矢量
r = sqrt(sum(sat2ref.^2, 2));               %矢量长度vector length
satind = ((sat2ref * [refpX; refpY; refpZ]) > 0) & (satdata.sathl == 0);%参考点处可见卫星satellites in view
% MEO_IGSO = (PRNNUM > 5) & (PRNNUM < 59);    %MEO/IGSO卫星
if satnum == 9
    satind([4, 5, 9, 14, 15, 19, 23, 30, 32, 34, 35]) = false;
elseif satnum == 12
    satind([14, 15, 19, 23, 30, 32, 34, 35]) = false;
elseif satnum == 15
    satind([23, 30, 32, 34, 35]) = false;
end
prnnum = PRNNUM(satind);                    %参考点处可见卫星编号satellites' prn in view
satlabels = mat2cell(satdata.prn(satind, :), ones(sum(satind), 1), 3);
%% 测量误差/检测统计量/门限
%测量误差variance of pseudorange measurement noise
sigmaURE = 4;
%检测统计量阈值the SSE threshold of diffetent alarm rate
alarmRate = [0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001];
SSE = zeros(length(alarmRate), 17);
alpha = 0:0.01:60;
for m = 1:length(alarmRate)
    for n = 1:17
        [B, I] = min(abs(cdf('Chisquare', alpha, n) - (1-alarmRate(m))));
        SSE(m, n) = alpha(I)*sigmaURE^2;
    end
end
%检测统计量门限判定，1超过门限
SPAlarm = @(x, n, m) x'*x > SSE(m, n);
%% 轨迹/卫星/真实距离
%位置
pos_r = [refpX; refpY; refpZ];%real position
%卫星坐标
satdata = BDSSatPosition(navdata, [2020 1 1 1 1 0], 'ECEF');    %仿真时刻的卫星数据
satpos = [satdata.x(satind), satdata.y(satind), satdata.z(satind)]';
%真实距离
r = mynorm(satpos - pos_r, 2, 1)';%satellites' geometric distance
%地球旋转修正add earth rotation error
earth_rot = 7.292115e-5 * (satpos(1, :) * pos_r(2) - satpos(2, :) * pos_r(1))' / c;
r = r + earth_rot;
%% MontoCarlo仿真
%欺骗卫星组合all combinations of satellites (the satellite number no less than 5)
spSatInd = logical(dec2bin(1:2^satnum-2, satnum) - '0');
spSatNum = sum(spSatInd, 2);
[spSatNum, ind] = sort(spSatNum, 'ascend');
spSatInd = spSatInd(ind, :);
spSatL2 = spSatNum <= spSatNmax;
spSatInd = spSatInd(spSatL2, :);
spSatNum = spSatNum(spSatL2);
% 进度条progress bar
tElapsed = 0;
spSatIdNum = length(spSatNum);
simTotNum = N*spSatIdNum;
simCurNum = 0;
simRemNum = simTotNum;
waitbarWin = waitbar(0, '初始化');
%分离判定the result of the algorithm
detect = zeros(spSatIdNum, N);       %分类结果是否被正确鉴定
detect0 = zeros(spSatIdNum, N);       %分类结果是否被正确鉴定
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
        % 进度条progress bar
        simCurNum = N*(spSatId-1)+n;
        simRemNum = simTotNum - simCurNum;
        waitbar(simCurNum/simTotNum, waitbarWin, ...
            ['仿真进度: ', num2str(simCurNum/simTotNum*100, '%.2f'), ' % (', ...
            num2str(simCurNum), '/', num2str(simTotNum), ')', ...
            '剩余时间:', num2str(remHour), ':', num2str(remMin, '%02d'), ':', num2str(remSec, '%02d')]);
        %% 伪距观测量生成
        %接收机时钟误差time bias
        dtr = randn()*1e-4;
        %接收伪距pseudorange
        pr = r;
        pr(spInd) = pr(spInd) + prangeBias(1) + (prangeBias(2)-prangeBias(1))*rand(spn, 1);
        URE = sigmaURE*randn(satnum, 1);
        pr = pr + c*dtr + URE;
        %定位解算pvt
        pos = [pos_r; 0]';
%         [pos, rtmp, pErr, G, H, cnt] = pvtpos0(pr(~spInd), satpos(:, ~spInd)', pos, 1e-10, 10);
        [pos, rtmp, pErr, G, H, cnt] = pvtpos0(pr, satpos', pos, 1e-10, 10);
        %% dual-RAIM
        spDet = false(satnum, 1);   %欺骗检测标记fault satellite flag
        spDet0 = false(satnum, 1);   %欺骗检测标记fault satellite flag
        if pErr'*pErr > SSE(alarmRateInd, satnum-4)
            % 检测到欺骗detect the existance of failure
            S = eye(satnum) - G*H*G';
            Fi = pErr.^2./diag(S);
            [Fimax, FiInd] = max(Fi);
            SSEi = pErr'*pErr - Fimax;
            if (SSEi < SSE(alarmRateInd, satnum-5))
                % 检测到仅存在单个异常伪距find only one fault
                spDet(FiInd) = true;
            else
                % 检测到多个异常伪距find more than one fault
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
                    % 检测到2个异常伪距find two faults
                    spDet([FijInd1, FijInd2]) = true;
                else
                    % 检测到多于2个异常伪距find more than two faults
                    spDet(:) = true;
                end
            end
            %% 精确计算剔除卫星后的SSE
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
        %% 鉴别情况分类
        falseNum(spSatId, n) = sum(xor(spDet, spInd));
        if all(~spDet)
            % 未检测到欺骗no faults
            detect(spSatId, n) = 0;
        elseif all(spDet == spInd)
            % 成功鉴别出所有欺骗find all faults
            detect(spSatId, n) = 1;
        elseif all(spDet)
            % 认为存在多于2颗异常卫星think there are more than two faults
            detect(spSatId, n) = 2;
        elseif any(spDet(~spInd)) && any(spInd(~spDet))
            % 将正常卫星误判为异常卫星并且有异常卫星没有检测出来
            % a right measurement detect as a fault and there is a fault
            % not detected
            detect(spSatId, n) = 3;
        elseif any(spDet(~spInd))
            % 将正常卫星误判为异常卫星
            % a right measurement detect as a fault
            detect(spSatId, n) = 4;
        elseif any(spInd(~spDet))
            % 有异常卫星没有检测出来
            % there is a fault not detected
            detect(spSatId, n) = 5;
        end
        %% 鉴别情况分类1
        falseNum0(spSatId, n) = sum(xor(spDet0, spInd));
        if all(~spDet0)
            % 未检测到欺骗
            detect0(spSatId, n) = 0;
        elseif all(spDet0 == spInd)
            % 成功鉴别出所有欺骗
            detect0(spSatId, n) = 1;
        elseif all(spDet0)
            % 认为存在多于2颗异常卫星
            detect0(spSatId, n) = 2;
        elseif any(spDet0(~spInd)) && any(spInd(~spDet0))
            % 将正常卫星误判为异常卫星并且有异常卫星没有检测出来
            detect0(spSatId, n) = 3;
        elseif any(spDet0(~spInd))
            % 将正常卫星误判为异常卫星
            detect0(spSatId, n) = 4;
        elseif any(spInd(~spDet0))
            % 有异常卫星没有检测出来
            detect0(spSatId, n) = 5;
        end
        n = n + 1;
    end
    tElapsed = toc;
end
close(waitbarWin);
toc;
%% 数据保存save data
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
%% 估算SSE分析
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
fprintf('平均成功概率：%.7f\t平均误检概率：%.7f\t平均失败概率：%.7f\n', 1-falseProb-failProb ,falseProb, failProb);

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
%% 精确计算SSE分析
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
fprintf('平均误检概率：%.4f\t平均失败概率：%.4f\n', falseProb, failProb);

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
