function [data_pad, NR] = pad_data(data, samplerate)
% pads a signal to reduce edge effects during filtering and resampling
% Inputs:
% data = nx1 column vector of data values (any units)
% samplerate = sample rate, in Hz
% Outputs:
% data_pad = mx1 column vector of padded data values (same units, same sample rate)
% NR = number of samples padded on either side of data

NR = samplerate;
before = 2*data(1) - flipud(data(2:NR+1)); % maintain continuous amplitude and slope
after = 2*data(end) - flipud(data(end-NR:end-1));
data_pad = [before; data; after]; % pad data

end