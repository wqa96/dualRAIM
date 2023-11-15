% ���㲻ͬ���������£�SSE�����龯���ʺ�©�����
clear;
% close all;

addpath data;
c = 299792458;              % ���٣�m/s��
%% ��������
N = 1000;
alarmRateInd = 2;
prangeBias = 100;
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
satind([14, 15, 19, 23, 30, 32, 34, 35]) = false;
prnnum = PRNNUM(satind);                    %�ο��㴦�ɼ����Ǳ��
satlabels = mat2cell(satdata.prn(satind, :), ones(sum(satind), 1), 3);
%% �ɼ���������ͼ
[satE, satN, satU] = ECEF2ENU1(sat2ref(satind, 1), sat2ref(satind, 2), sat2ref(satind, 3), ...
    refpLLA(1), refpLLA(2));                %�ɼ����ǵ�ENU����
sat0 = [satE, satN, satU];
phi = asin(satU ./ r(satind)) * 180 / pi;   %�ɼ����ǵ�����
theta = atan2(satE, satN) * 180 / pi;       %�ɼ����ǵķ�λ��
r = sqrt(satE.^2 + satN.^2 + satU.^2);      %�ɼ����ǵ���ʵ����
% figure(11);
% skyPlot0(theta, phi, prnnum);
%% α��/���ͳ����
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
%��λ���
pos_r = zeros(1, 3);
pos_s0 = zeros(1, 3);
%��������
satdata = BDSSatPosition(navdata, t, 'ECEF');   %����ʱ�̵���������
sat2ref = [satdata.x(satind) - refpX, ...
    satdata.y(satind) - refpY, ...
    satdata.z(satind) - refpZ];                     %����λ�ù��ڲο�λ�õ�ʸ��
[satE, satN, satU] = ECEF2ENU1(sat2ref(:, 1), sat2ref(:, 2), sat2ref(:, 3), ...
    refpLLA(1), refpLLA(2));        %�ɼ����ǵ�ENU����
sat = [satE, satN, satU];
%��ʵα��
r = sqrt(sum((sat - pos_r).^2, 2));
%% Monte Carlo ����
%���д��ڵ���7���ǵ����
satnum = 12;
satInd = logical(dec2bin(1:2^satnum-1)-'0');
satnum = sum(satInd, 2);
satInd = satInd(satnum > 5, :);
satnum = satnum(satnum > 5);
Num = length(satnum);
% �۲���
missDet1 = nan(Num, 12);        %���쳣SSE���©��
missDet1all = nan(Num, 12, N);
missDet2 = nan(Num, 66);        %˫�쳣SSE���©��
missDet2all = nan(Num, 66, N);
Sii = nan(Num, 12);             %���쳣��������
Fi = nan(Num, 12, N);
Pij = nan(Num, 66, 2);          %˫�쳣
Sij = nan(Num, 66, 3);
Fij = nan(Num, 66, N);
% ������
tElapsed = 0;
simTotNum = Num;
simCurNum = 0;
simRemNum = simTotNum;
waitbarWin = waitbar(0, '��ʼ��');
tic;
for num = 1:Num
    nn = 1;
    % ������
    simCurNum = num-1;
    simRemNum = simTotNum - simCurNum;
    simSpeed = simCurNum/tElapsed;
    remSec = simRemNum/simSpeed;
    remHour = floor(remSec/3600);
    remMin = floor(mod(remSec, 3600)/60);
    remSec = floor(mod(remSec, 60));
    waitbar(simCurNum/simTotNum, waitbarWin, ...
        ['�������: ', num2str(simCurNum/simTotNum*100, '%.2f'), ' % (', ...
        num2str(simCurNum), '/', num2str(simTotNum), ')', ...
        'ʣ��ʱ��:', num2str(remHour), ':', num2str(remMin, '%02d'), ':', num2str(remSec, '%02d')]);
    % ��ǰ��������ǲ���
    r1 = r(satInd(num, :));
    sat1 = sat(satInd(num, :), :);
    %���ջ�ʱ�����
    dtr = randn()*1e-7;
    pr = r1 + c*dtr;
    %% �������е��쳣�������
    for kk = 1:satnum(num)
        % �������ǵ�GDOP
        flag = true(satnum(num), 1); flag(kk) = false;
        sattmp = sat1(flag, :);
        rtmp = r1(flag);
        G = -(sattmp - pos_r)./rtmp;
        G(:, 4) = 1;
        H = inv(G'*G);
        g1 = [-(sat1(~flag, :) - pos_r)./r1(~flag), 1]';
        Sii(num, kk) = 1/(1 + g1'*H*g1);
        % ���ƫ����
        pr1 = pr;
        pr1(kk) = pr1(kk) + prangeBias;
        misscnt = 0;
        %% N�η���
        for nn = 1:N
            % �������
            URE = sigmaURE*randn(satnum(num), 1);
            pr1nn = pr1 + URE;
            % ��λ����
            ptmp = [pos_r, 0];
            [pos, rtmp, perr, G, H, cnt] = positioning(pr1nn, sat1, ptmp, 1e-10, 10);
            sseFlag = perr'*perr < SSE(alarmRateInd, satnum(num)-4);
            if sseFlag
                misscnt = misscnt + 1;
                missDet1all(num, kk, nn) = 1;
            else
                missDet1all(num, kk, nn) = 0;
            end
            Fi(num, kk, nn) = Sii(num, kk)*(prangeBias+URE(~flag))^2;
        end
        missDet1(num, kk) = misscnt/N;
    end
    %% ��������˫�쳣�������
    if satnum(num) < 7
        continue;
    end
    flags = nchoosek(1:satnum(num), 2);
    for kk = 1:size(flags, 1)
        % �������ǵ�GDOP
        flag = true(satnum(num), 1); flag(flags(kk, :)) = false;
        sattmp = sat1(flag, :);
        rtmp = r1(flag);
        G = -(sattmp - pos_r)./rtmp;
        G(:, 4) = 1;
        H = inv(G'*G);
        g1 = [-(sat1(~flag, :) - pos_r)./r1(~flag), [1; 1]]';
        P = eye(2) + g1'*H*g1;
        [V, D] = eig(P);
        Pij(num, kk, :) = diag(D);
        Sij(num, kk, 1:2) = diag(P);
        Sij(num, kk, 3) = P(1, 2);
        % ���ƫ����
        pr1 = pr;
        pr1(~flag) = pr1(~flag) + prangeBias;
        misscnt = 0;
        %% N�η���
        for nn = 1:N
            % �������
            URE = sigmaURE*randn(satnum(num), 1);
            pr1nn = pr1 + URE;
            % ��λ����
            ptmp = [pos_r, 0];
            [pos, rtmp, perr, G, H, cnt] = positioning(pr1nn, sat1, ptmp, 1e-10, 10);
            sseFlag = perr'*perr < SSE(alarmRateInd, satnum(num)-4);
            if sseFlag
                misscnt = misscnt + 1;
            end
        end
        missDet2(num, kk) = misscnt/N;
    end
    tElapsed = toc;
end
close(waitbarWin);
%% ���ݱ���
t = datetime('now');
tMDHM = reshape(num2str([t.Month, t.Day, t.Hour, t.Minute]', '%02d')', 1, []);
fullpath = mfilename('fullpath');
[~, name] = fileparts(fullpath);
simMode = num2str(prangeBias');
savefilename = ['.\data\', name, '_', tMDHM, '_', simMode, '_', num2str(alarmRateInd), '.mat'];
save(savefilename);
%% ���쳣���Ƿ���
satNumMin = 6;
% load('SSETestVsSatnum_03301954_100_1.mat');
figid = 1;
% figure(figid); figid = figid + 1;
% tmp = mean(missDet1, 2, 'omitnan');
% plot(satnum, tmp, 'o');
% xlabel('��������'); ylabel('�龯����');

% figure(figid); figid = figid + 1;
% tmp = mean(missDet1, 2, 'omitnan');
% missRate = zeros(12-satNumMin+1, 1);
% for n = satNumMin:12
%     missRate(n-satNumMin+1) = mean(tmp(satnum == n));
% end
% semilogy(satNumMin:12, missRate); xlim([satNumMin, 12]);
% xlabel('��������'); ylabel('©�����');

% figure(figid); figid = figid + 1;
% tmp = mean(Sii, 2, 'omitnan');
% plot(satnum, tmp, 'o');
% xlabel('��������'); ylabel('S_{i,i}');

% figure(figid); figid = figid + 1;
% tmp = mean(Sii, 2, 'omitnan');
% sii = zeros(12-satNumMin+1, 1);
% for n = satNumMin:12
%     sii(n-satNumMin+1) = mean(tmp(satnum == n));
% end
% plot(satNumMin:12, sii); xlim([satNumMin, 12]);
% xlabel('��������'); ylabel('S_{i,i}');

figure(figid); figid = figid + 1;
loglog(Sii(:), missDet1(:), '.');
xlabel('S_{i,i}'); ylabel('©�����');

figure(figid); figid = figid + 1;
for nn = 6:12
    flag = satnum == nn;
    xtmp = Sii(flag, :);
    ytmp = missDet1(flag, :);
    loglog(xtmp(:), ytmp(:), '.'); hold on;
end
hold off;
xlabel('S_{i,i}'); ylabel('©�����');
%% ˫�쳣���Ƿ���
% figure(figid); figid = figid + 1;
% tmp = mean(missDet2, 2, 'omitnan');
% missRate = zeros(12-satNumMin+1, 1);
% for n = satNumMin:12
%     missRate(n-satNumMin+1) = mean(tmp(satnum == n));
% end
% semilogy(satNumMin:12, missRate); xlim([satNumMin, 12]);
% xlabel('��������'); ylabel('©�����');

% tmp1 = Pij(:, :, 1); tmp2 = Pij(:, :, 2);
% figure(figid); figid = figid + 1;
% plot3(tmp1(:), tmp2(:), missDet2(:), '.');
% grid on;
% set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log');

% figure(figid); figid = figid + 1;
% plot3(1./tmp1(:), 1./tmp2(:), missDet2(:), '.');
% grid on;
% set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log');

% figure(figid); figid = figid + 1;
% loglog(tmp1(:)+tmp2(:), missDet2(:), '.');

% figure(figid); figid = figid + 1;
% loglog(1./tmp1(:)+1./tmp2(:), missDet2(:), '.');

% tmp1 = Sij(:, :, 1); tmp2 = Sij(:, :, 2); tmp3 = Sij(:, :, 3);
% figure(figid); figid = figid + 1;
% plot3(tmp1(:), tmp2(:), missDet2(:), '.');
% grid on;
% set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log');

% figure(figid); figid = figid + 1;
% plot3(tmp1(:), tmp2(:), tmp3(:), '.');
% grid on;

% flag = abs(tmp3) < 100;
% figure(figid); figid = figid + 1;
% plot3(tmp1(flag), tmp2(flag), missDet2(flag), '.');
% grid on;
% set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log');

figure(figid); figid = figid + 1;
tmp = (Sij(:, :, 1) + Sij(:, :, 2) - 2*Sij(:, :, 3))./(Sij(:, :, 1).*Sij(:, :, 2)-Sij(:, :, 3).^2);
loglog(tmp(:), missDet2(:), '.'); grid on;
%% �龯����Vs©�����
% satnum = 12;
% satInd = logical(dec2bin(1:2^satnum-1)-'0');
% satnum = sum(satInd, 2);
% satInd = satInd(satnum > 4, :);
% satnum = satnum(satnum > 4);
% Num = length(satnum);
% 
% % datafiles = {'SSETestVsSatnum_03301819_100_1.mat';
% %     'SSETestVsSatnum_03301843_100_2.mat';
% %     'SSETestVsSatnum_03301847_100_3.mat';
% %     'SSETestVsSatnum_03301851_100_4.mat';
% %     'SSETestVsSatnum_03301855_100_5.mat';
% %     'SSETestVsSatnum_03301859_100_6.mat';
% %     };
% datafiles = {'SSETestVsSatnum_03301954_100_1.mat';
%     'SSETestVsSatnum_03301954_100_2.mat';
%     'SSETestVsSatnum_03301954_100_3.mat';
%     'SSETestVsSatnum_03301955_100_4.mat';
%     'SSETestVsSatnum_03301954_100_5.mat';
%     'SSETestVsSatnum_03301955_100_6.mat';
%     };
% missDets = cell(6, 1);
% for k = 1:6
%     tmp = load(datafiles{k}, 'missDet');
%     missDets{k} = tmp.missDet;
% end
% 
% missRates = zeros(7, 6);
% for k = 1:6
%     tmp = mean(missDets{k}, 2, 'omitnan');
%     for n = satNumMin:12
%         missRates(n-satNumMin+1, k) = mean(tmp(satnum == n));
%     end
% end
% 
% figid = 1;
% figure(figid); figid = figid + 1;
% semilogy(satNumMin:12, missRates); xlim([satNumMin, 12]);
% xlabel('��������'); ylabel('©�����');
% legend({'p_{alarm}=0.1';
%     'p_{alarm}=0.05';
%     'p_{alarm}=0.01';
%     'p_{alarm}=0.001';
%     'p_{alarm}=0.0001';
%     'p_{alarm}=0.00001';});
% 
% figure(figid); figid = figid + 1;
% semilogx([0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001], ...
%     missRates(1, :)*100);
% xlabel('�龯����'); ylabel('©�����(%)');
%% �龯����Ϊ0.05��©�������Sii�Ĺ�ϵ
alarmRate = 0.05;
TH = zeros(1, 17); %��ͬ���ɶȵ���ֵ
alpha = 0:0.01:60;
for n = 1:17
    missRates = cdf('Chisquare', alpha, n);
    [B, I] = min(abs((1-missRates) - alarmRate));
    TH(n) = alpha(I);
end

Bi = 100;
Sii0 = (0:1e-3:1)';
missRates1 = zeros(length(Sii0), 17);
for n = 1:16
    dx = 0.01;
    x = -50:dx:150;
    th = (16*TH(n+1) - x.^2)/16;
    px = pdf('Normal', x, Bi*sqrt(Sii0), 4);
    missRates1(:, n) = sum(cdf('Chisquare', th, n) .* px * dx, 2);
end

for nn = 6:6
    flag = satnum == nn;
    xtmp = Sii(flag, :);
    ytmp = missDet1(flag, :);
    figure; loglog(1./xtmp(:), ytmp(:), '.');
    hold on; loglog(1./Sii0, missRates1(:, nn-5)); hold off;
    ylim([1e-3, 1]);
    xlabel('$1+g_i^TH_ig_i$', 'Interpreter', 'latex');
    ylabel('Miss detection rate');
    legend('simulation result', 'ideal');
end
%% �龯����Ϊ0.05��©�������Qij�Ĺ�ϵ
alarmRate = 0.05;
TH = zeros(1, 17); %��ͬ���ɶȵ���ֵ
alpha = 0:0.01:60;
for n = 1:17
    missRates = cdf('Chisquare', alpha, n);
    [B, I] = min(abs((1-missRates) - alarmRate));
    TH(n) = alpha(I);
end
b = 100;
lambda0 = (0:1e-3:1)';
sigma = 4;
missRates2 = zeros(length(lambda0), 17);
for n = 1:15
    lambda = lambda0*b^2/sigma^2;
    dx = 0.0001;
%     x2 = ones(length(lambda), 1)*(0:dx:25*25*2);
    meanBij2 = lambda + 2;
    varBij2 = 4*lambda + 4;
    x2 = meanBij2 - 10*sqrt(varBij2) + 20*sqrt(varBij2)*(0:dx:1);
    th = TH(n+2) - x2;
%     px = ncx2pdf(x2, 2, lambda*ones(1, size(x2, 2)));
    px = ncx2pdf(x2, 2, lambda*ones(1, size(x2, 2)));
    delta = 20*sqrt(varBij2) * dx;
    missRates2(:, n) = sum(cdf('Chisquare', th, n) .* px .* delta, 2);
end
for nn = 7:12
    flag = satnum == nn;
    xtmp = (Sij(flag, :, 1) + Sij(flag, :, 2) - 2*Sij(flag, :, 3))./(Sij(flag, :, 1).*Sij(flag, :, 2)-Sij(flag, :, 3).^2);
    ytmp = missDet2(flag, :);
    figure; loglog(xtmp(:), ytmp(:), '.'); grid on;
    hold on; loglog(lambda0, missRates2(:, nn-6)); hold off;
    ylim([1e-3, 1]);
    xlabel('$\lambda_0$', 'Interpreter', 'latex');
    ylabel('Miss detection rate');
    legend('simulation result', 'ideal');
end

