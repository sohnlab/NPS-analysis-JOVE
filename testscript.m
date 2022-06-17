%% testing

clearvars; close all;

filepath_JOVE12 = fullfile('example data', '20220617_A549_dev2B_w12_p25_try1.mat');
filepath_JOVE10 = fullfile('example data', '20220616_BEAS2B_Dev1B_w10_p30_try1.mat');
filename_04 = fullfile('example data', 'EDL_11_04_19_test1_2018.mat');
filename_21 = fullfile('example data', 'EDL_P1_w_11_21_test2_cropped.mat');
filename_fs = fullfile('example data', 'sample_rate.mat');

%% using chondrocyte example data

% % sample rate: assumed to be 50kHz
% load(filename_fs, 'fs'); % [Hz]
% 
% % effective diameter
% %   > here, just using the value from the original version of
% %       mNPS_readJOVE.m (March 2020), which also set h=18.21
% %   > for accurate analysis, this value needs to be determined by a calibration
% %       experiment for the wafer that was used (depends on mask design & wafer height)
% %   > D_e for the contraction segment is calculated in mNPS_readJOVE.m based on
% %       De_np and the wafer design
% De_np = 22.5; % D_e for reference & recovery segments [um]
% 
% % 11/04/2019 (chondrocytes)
% ch_height = 18.21; % [um] for devices used before 11/15/2019
% thresholds = [3e-4, 1e-3];
% out_04 = mNPS_procJOVE(filename_04, ch_height, De_np, thresholds, fs);
% 
% % % 11/21/2019 test2 (P1 chondrocytes w/ GF's)
% % ch_height = 16.6; % [um] for devices used on or after 11/15/2019
% % out_21 = mNPS_procJOVE(filename_21, ch_height, De_np, [], fs);

%% using JOVE example data

% sample rate: saved in file
S_in = load(filepath_JOVE12);
fs = S_in.sampleRate;

% effective diameter: estimate based on geometry
%   > for accurate analysis, this value needs to be determined by a calibration
%       experiment for the wafer that was used (depends on mask design & wafer height)
%   > D_e for the contraction segment is calculated in mNPS_readJOVE.m based on
%       De_np and the wafer design
w_NP = 25; % channel width in node-pore sections [um]

% 06/17/2022 (A549, wC=12)
ch_height = 30; % [um] (approximate)
De_np = 2 * ch_height * w_NP / (ch_height + w_NP); % effective diameter estimate
thresholds = []; % use defaults
out_04 = mNPS_procJOVE(filepath_JOVE12, ch_height, De_np, thresholds, fs);

% 06/16/2022 (BEAS2B, wC=10)
ch_height = 30; % [um] (approximate)
De_np = 2 * ch_height * w_NP / (ch_height + w_NP); % effective diameter estimate
thresholds = []; % use defaults
out_04 = mNPS_procJOVE(filepath_JOVE10, ch_height, De_np, thresholds, fs);
