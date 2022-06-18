%% testing

clearvars; close all;

filepath_JOVE12 = fullfile('example data', '20220617_A549_dev2B_w12_p25_try1.mat');
filepath_JOVE10 = fullfile('example data', '20220616_BEAS2B_Dev1B_w10_p30_try1.mat');

%% parameters for JOVE example data

% effective diameter: estimate based on geometry
%   > for accurate analysis, this value needs to be determined by a calibration
%       experiment for the wafer that was used (depends on mask design & wafer height)
%   > D_e for the contraction segment is calculated in mNPS_readJOVE.m based on
%       De_np and the wafer design
w_NP = 25; % channel width in node-pore sections [um]

%% 06/17/2022 (A549, wC=12)

% sample rate: saved in file
S_in = load(filepath_JOVE12, 'sampleRate');
fs = S_in.sampleRate;

wC = 12; % [um] contraction channel width
ch_height = 30; % [um] (approximate)
De_np = 2 * ch_height * w_NP / (ch_height + w_NP); % effective diameter estimate

ASLS_param = struct('lambda',1e10, 'p',0, 'noise_margin',2.4e-4, 'max_iter',40);
thresholds = [1e-4, 5e-4];

out_j12 = mNPS_procJOVE(filepath_JOVE12, ch_height, De_np, wC, thresholds, fs, ASLS_param);

%% 06/16/2022 (BEAS2B, wC=10)

% % sample rate: saved in file
% S_in = load(filepath_JOVE10, 'sampleRate');
% fs = S_in.sampleRate;
% 
% wC = 10; % [um] contraction channel width
% ch_height = 30; % [um] (approximate)
% De_np = 2 * ch_height * w_NP / (ch_height + w_NP); % effective diameter estimate
% 
% ASLS_param = struct('lambda',1e10, 'p',0, 'noise_margin',1.8e-4, 'max_iter',40);
% thresholds = [0.9e-4, 4e-4];
% 
% out_j10 = mNPS_procJOVE(filepath_JOVE10, ch_height, De_np, wC, thresholds, fs, ASLS_param);
