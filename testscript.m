%% Kristen testing

clearvars;
filename_04 = 'example data/EDL_11_04_19_test1_2018.mat';
filename_21 = 'example data/EDL_P1_w_11_21_test2_cropped.mat';
filename_fs = 'example data/sample_rate.mat';

% sample rate: assumed to be 50kHz
load(filename_fs, 'fs'); % [Hz]

% effective diameter
% here, just using the value from the original version of
%   mNPS_readChon.m (March 2020), which also set h=18.21
% for accurate analysis, this value needs to be determined by a calibration
%   experiment for the wafer that was used (depends on mask design & wafer height)
% D_e for the contraction segment is calculated in mNPS_readChon.m based on
%   De_np and the wafer design
De_np = 22.5; % D_e for reference & recovery segments [um]

%% test mNPS_procChon

% 11/04/2019
ch_height = 18.21; % [um] for devices used before 11/15/2019
thresholds = [3e-4, 1e-3];
out_04 = mNPS_procChon(filename_04, ch_height, De_np, thresholds, fs);

% % 11/21/2019 test2 (P1 chondrocytes w/ GF's)
% ch_height = 16.6; % [um] for devices used on or after 11/15/2019
% out_21 = mNPS_procChon(filename_21, ch_height, De_np, [], fs);
