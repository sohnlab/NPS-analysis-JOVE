function data = unpad_data(data_pad, samplerate)
% undoes the effects of pad_data

data = data_pad(samplerate+1:end-samplerate);

end