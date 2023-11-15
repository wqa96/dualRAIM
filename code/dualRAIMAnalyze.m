clear;
close all;

datafiles = {'dualRAIMSim_1202_1_10000m_04031635.mat';
    'dualRAIMSim_1202_2_10000m_04031635.mat';
    'dualRAIMSim_1202_3_10000m_04031635.mat';
    'dualRAIMSim_1202_4_10000m_04031635.mat';
    'dualRAIMSim_1202_5_10000m_04031635.mat';
    'dualRAIMSim_1202_6_10000m_04031635.mat'};

for n = 1:6
    data(n) = load(datafiles{n});
end

