function [y_smoothed, y_downsampled, y_detrended] = mNPS_fastQC(data, sampleRate, k_sample, detrend_flag, filt_flag, ASLS_param)
% Fast option for post-processing mNPS signal.
%   default sampleRate = 50000 [Hz]
%   default k_sample (downsample factor) = 20
%   optional bool detrend_flag to perform baseline subtraction (default=false)
%   optional bool filt_flag to enable 60 Hz bandstop filter (default=false)
%   optional struct ASLS_param (baseline fitting parameters passed to ASLS.m)

    %% parse inputs
    
    % if data is a 2xn matrix, just extract the second row (voltage)
    if size(data,1)==2 && size(data,2)>2
        data = data(2,:)';
    % if data is a row vector, transpose it to a column vector
    elseif size(data,1)==1 && size(data,2)>1
        data = data';
    % the only other acceptable form for data is a column vector
    elseif ~( size(data,2)==1 && size(data,1)>1 )
        error('data must be a vector or a 2xn array');
    end
    
    if nargin<5 || isempty(filt_flag)
        filt_flag = false;
    end

    if nargin<4 || isempty(detrend_flag)
        detrend_flag = false;
    end

    if nargin<3 || isempty(k_sample)
        k_sample = 20; % to downsample from 50 kHz to 2.5 kHz
    end

    if nargin<2 || isempty(sampleRate)
        sampleRate = 50e3; % 50 kHz
    end

    % default ASLS parameters (if needed)
    if ( nargin<6 || isempty(ASLS_param) ) && detrend_flag
        ASLS_param = struct();
        ASLS_param.lambda = 1e9; % default 1e5; larger=smoother, smaller=wiggly-er (may not be unit-independent)
        ASLS_param.p = 3e-3; % default 0.01; 0>p>1 (as low as possible while still converging)
        ASLS_param.noise_margin = 1e-4; % default 0; allows baseline to sit within the baseline noise
        ASLS_param.max_iter = 20; % default 5; just make sure it converges
    end
    
    %% perform preprocessing
    
    % invert data if needed
    if (mean(data) < 0 )
        data = -data;
    end
    
    % perform rectangular smoothing
    y_smoothed = fastsmooth(data,200,1,1);
    
    % apply 60Hz badstop filter if desired
    if filt_flag
        d = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
            'DesignMethod','butter','SampleRate',sampleRate);
        y_smoothed = filtfilt(d,y_smoothed);
    end
    
    % downsample by period N
    y_downsampled = downsample(y_smoothed,k_sample);
    
    % subtract baseline if desired
    if detrend_flag
        y_baseline = -1 * ASLS(-1*y_downsampled, ASLS_param);
        y_detrended = y_downsampled - y_baseline;
    else
        y_detrended = [];
    end
    
    %% plot results
    
    figure();
    
    if detrend_flag
        plot(y_detrended);
        titlestr = 'data: smoothed, downsampled, and detrended';
    else
        plot(y_downsampled);
        titlestr = 'data: smoothed & downsampled';
    end
    
    if filt_flag
        title([titlestr,' (with 60Hz bandstop filter)']);
    else
        title(titlestr);
    end
    
end

