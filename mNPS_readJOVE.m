function [OUT_array, empty, auto_thresh_value, column_names, column_units, rec_cat_description] = ...
    mNPS_readJOVE(data_vector, sampleRate, ch_height, De_np, wC, thresholds, plotflag, fitflag, ASLS_param)
    % Reads mNPS data and returns OUT_array matrix
    % INPUTS:
    %   data_vector = row vector of doubles
    %   sampleRate = int [Hz]
    %   ch_height = double [um]: channel height measured from SU-8 wafer
    %   De_np = double [um]: effective diameter for the node-pore segments (from REF device)
    %   wC = double [um]: width of the contraction segment
    %   thresholds = 1x2 vector of doubles: [low_threshold, high_threshold]
    %   plotflag = bool: whether to plot the window data
    %   fitflat = bool: whether to perform mNPS-r recovery curve fitting
    %   ASLS_param = struct (optional): baseline fitting parameters passed to ASLS.m

    if nargin<9
        ASLS_param = []; % use defaults
    end

    %% SECTION 0: device parameters

    % segment layout
    total_segs = 7;
    num_ref_segs = 3;
    num_rec_segs = 3;
    % check layout
    if total_segs ~= (num_ref_segs + num_rec_segs + 1)
        error('check device geometry! number of segments doesn''t match!');
    end

    % mask geometry [um]
    L = 8130; % sNPS_ver2.1 (both) - total length of the NPS channel (from start of 1st pore to end of last pore)
    npL_ref = [1155, 1155, 577.5]; % sNPS_ver2.1 (both)
    npL_rec = [577.5, 1155, 1155]; % sNPS_ver2.1 (both)
    sqL = 2055; % sNPS_ver2.1 (both)
    wNP = 25; % sNPS_ver2.1 (both)

    % calculate De_c based on De_np
    De_c = De_np*(wC/wNP)^(0.5); % D_e (effective diameter) for contraction segment

    %% SECTION 1: load data and perform basic signal conditioning

    Fs = sampleRate/1000; % convert to kHz

    N = 20; % downsample factor

    % flip if negative
    if (sum(data_vector<0) / length(data_vector)) > 0.5
        data_vector = -data_vector;
        warning('inverting raw data because the majority of datapoints were negative');
    end

    y_smoothed = fastsmooth(data_vector',200,1,1); % perform rectangular smoothing

    ym = downsample(y_smoothed,N); % downsample by period N

    if size(ym,1) < size(ym,2)
        ym = ym'; % transpose if vector is of the wrong dimension
    end

    % remove baseline
    y_baseline = -1 * ASLS(-1*ym, ASLS_param);
    y_detrend = ym - y_baseline;

    %% SECTION 2: threshold signal by differences
    % take the difference of ym, threshold by lower value
    % thresholds for differences, user provided

    ym_diff = diff(y_detrend); % compute difference
    ym_diff(abs(ym_diff) < thresholds(1)) = 0; % threshold values below thresholds(1)
    ym_diff(1:10) = 0; % zero out the first few values
    ym_diff(end-10:end) = 0; % zero out the last few values

    % ensure all values in squeeze channel are zero
    squeeze_begin = find(ym_diff <= -thresholds(2),1); % find where squeeze starts

    % advance by 1 until large difference ends
    while (ym_diff(squeeze_begin) <= -thresholds(2))
        squeeze_begin = squeeze_begin + 1;
    end

    squeeze_end = find(ym_diff(squeeze_begin+20:end) >= thresholds(2),1); % find end of squeeze
    ym_diff(squeeze_begin+1:squeeze_begin+squeeze_end+20) = 0; % zero out the squeeze channel

    %% SECTION 3: identify nonzero differences (nz_mat is the matrix of nonzero differences)

    nz_mat = ones(2,length(nonzeros(ym_diff))); % preallocation
    k = 1; % index for A
    for i = 1:length(ym_diff)
        if (ym_diff(i) ~= 0) % look for nonzeros
            nz_mat(1,k) = i; % array index of nonzero
            nz_mat(2,k) = ym_diff(i); % nonzero value
            k = k+1;
        end
    end

    %% SECTION 4: remove error from A

    k=1;
    while (k < ceil(log(length(nz_mat))))
        i=1;
        while (i < length(nz_mat))

            % Case 1: current and next both positive && next > current
            if nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                    nz_mat(2,i+1) > nz_mat(2,i)
                % move next into current
                nz_mat(:,i) = nz_mat(:,i+1);

            % Case 2: current and next both positive && current > next
            elseif nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                    nz_mat(2,i) > nz_mat(2,i+1)
                % move current into next
                nz_mat(:,i+1) = nz_mat(:,i);

            % Case 3: current and next both negative && current > next
            elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                    nz_mat(2,i) > nz_mat(2,i+1)
                % move next into current
                nz_mat(:,i) = nz_mat(:,i+1);

            % Case 4: current and next both negative && next > current
            elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                    nz_mat(2,i+1) > nz_mat(2,i)
                % move current into next
                nz_mat(:,i+1) = nz_mat(:,i);
            end

            i=i+1;
        end
        k=k+1;
    end

    %% SECTION 5: remove repeats in A

    unique_xs = unique(nz_mat(1,:)); % unique values
    unique_is = ones(1,length(unique_xs)); % unique indices
    for i = 1:length(unique_is)
        unique_is(i) = find(nz_mat(1,:) == unique_xs(i), 1);
    end
    unique_ys = nz_mat(2,unique_is);
    nz_mat = [unique_xs; unique_ys];

    %% SECTION 6: rectangularize pulses

    ym_rect = y_detrend;
    k = 1;
    while (k <= 50)
        i = 1;
        while (i < length(nz_mat))

            if nz_mat(2,i) < 0 && nz_mat(2,i+1) > 0 % look for sign change in differences
                % replace all values in between with mean
                ym_rect(nz_mat(1,i):nz_mat(1,i+1)) = ...
                    mean(y_detrend(nz_mat(1,i):nz_mat(1,i+1))); 
            end

            i=i+1;
        end
        k=k+1;
    end

    %% SECTION 7: Plot figures if flag is true

    if max(ym_diff) < thresholds(1) % waste of time, no pulse
        plotflag = false;
        empty = true;
    else
        empty = false;
    end

    if plotflag == true
        Pix_SS = get(0,'screensize');
        figh = figure(42); 
        figsize = [0.1 0.1 0.45 0.75]*Pix_SS(4);
        set(figh,'units','pixels','pos',figsize);

        % take top and bottom 3 values
        nsorted_d = sort(ym_diff);
        min_vals = nsorted_d(1:3);

        psorted_d = sort(ym_diff,'descend');
        max_vals = psorted_d(1:3);

        % set auto-thresholds
        if abs(min_vals(3)) < abs(max_vals(3))
            auto_thresh_value = abs(min_vals(3));
        else
            auto_thresh_value = abs(max_vals(3));
        end

        % difference plot
        ax1 = subplot(3,1,1);
        figwin_tighten();
        difp = plot(ym_diff,'k-'); difp.LineWidth = 1;
        title('y_{diff}');
        set(gca,'FontSize',10);
        grid(ax1,'on');
        ax1.XMinorGrid = 'on';
        axis([0, length(ym_diff), 1.1*mean(min_vals), 1.1*mean(max_vals)]);
        for i = 1:length(max_vals)
            label_str = sprintf('%3.3e',min_vals(i));
            text(i*600,1.35*max_vals(1),label_str,'FontSize',12);

            label_str = sprintf('%3.3e',max_vals(i));
            text(i*600,1.85*max_vals(1),label_str,'FontSize',12);
        end

        % plot thresholds
        hold on
        linel = length(ym_diff);
        bthlin_l = line([0 linel], [-thresholds(1), -thresholds(1)]);
        bthlin_u = line([0 linel], [thresholds(1), thresholds(1)]);

        tthlin_l = line([0 linel], [-thresholds(2), -thresholds(2)]);
        tthlin_u = line([0 linel], [thresholds(2), thresholds(2)]);

        bthlin_l.Color = [1 0 0];
        bthlin_u.Color = [1 0 0];
        tthlin_l.Color = [0 0 1];
        tthlin_u.Color = [0 0 1];
        hold off

        % rectangularized
        ax2 = subplot(3,1,2);
        figwin_tighten();
        plot(ym_rect);
        title('y_{rect}');
        set(gca,'FontSize',10);
        grid(ax2,'on');
        ax2.XMinorGrid = 'on';
        axis([0, length(ym_rect), 1.1*min(ym_rect), 0.01]);

        % smoothed
        subplot(3,1,3);
        figwin_tighten();
        plot(y_smoothed,'k-');
        title('y_{LP}');
        set(gca,'FontSize',10);
        axis([0, length(y_smoothed), 0.999*min(y_smoothed), 1.001*max(y_smoothed)]);

    else
        auto_thresh_value = [];
    end

    %% SECTION 8: Detect NPS pulses
    % pulse_series is a matrix with the indices and parameters for rectangular pulses

    i=1;
    k = 0;
    backset = 10;
    pulse_series = ones(length(nz_mat),4);
    while (i < length(nz_mat))
        if nz_mat(2,i) < 0 && nz_mat(2,i+1) > 0 % starts negative and flips sign
            k = k + 1;
            pulse_series(k,1) = nz_mat(1,i); % Start index
            pulse_series(k,2) = nz_mat(1,i+1); % End index
            pulse_series(k,3) = mean(ym((nz_mat(1,i)-backset):(nz_mat(1,i)-backset+10))); % normalized baseline current
            pulse_series(k,4) = mean(y_detrend(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % avg current drop between pulses
            pulse_series(k,5) = std(y_detrend(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % std dev of current drop
        end
        i=i+1;
    end

    % remove empty entries in P
    cci = 1;
    stopc = size(pulse_series,1);
    while(cci <= stopc)
        if (pulse_series(cci,1) == 1 && pulse_series(cci,2) == 1)
            pulse_series(cci,:) = [];
            stopc = stopc - 1;
        else
            cci = cci + 1;
        end
    end

    %% SECTION 9: Extract mNPS pulse data

    % preallocate array of NPS pulse data
    out = nan(length(pulse_series) - (total_segs-1), 12);

    for k = 1 : length(pulse_series)+1-total_segs

        % get relevant row numbers
        sq_k = k + num_ref_segs; % contraction segment
        ref_k_start = k; % first reference segment
        ref_k_end = sq_k - 1; % last reference segment
        rec_k_start = sq_k + 1; % first recovery segment
        rec_k_end = rec_k_start + num_rec_segs - 1; % last recovery segment
        % check segment indices
        if rec_k_end ~= (k + total_segs - 1)
            error('check segment indexing!');
        end

        start_index = pulse_series(k,1); % starting index
        I_baseline = pulse_series(k,3); % baseline current
        
        % average dI & dT in reference segments

        % average node-pore current drop in reference segments
        dI_np = -mean(pulse_series(ref_k_start:ref_k_end,4));

        % average node-pore transit time in reference segments [ms]
        %   *** not applicable in JOVE device designs (sNPS_v2.1) bc the segments are of unequal lengths
        dT_np = nan;

        % node-pore transit time in each reference segment [ms] (row vector)
        dT_np_segs = ( pulse_series(ref_k_start:ref_k_end,2) - pulse_series(ref_k_start:ref_k_end,1) )' ./Fs.*N;

        % dI & dT in contraction (sqeeze) segment
        dI_c = -pulse_series(sq_k,4); % squeeze current drop
        dI_c_std = pulse_series(sq_k,5); % std. dev. of squeeze current drop
        dT_c = (pulse_series(sq_k,2) - pulse_series(sq_k,1)) /Fs*N; % squeeze transit time (ms)

        % post-squeeze NP current drops
        dI_rec = -pulse_series(rec_k_start:rec_k_end, 4);
        
        %% determine recovery time & category
        % recovery time is determined when post-squeeze NP current drop
        %   reaches pre-squeeze NP current drop (within 8% error threshold)
        % if the cell has "instant" recovery, recovery time is defined as
        %   elapsed time between the end of the squeeze segment and the
        %   beginning of the first recovery segment
        % if the cell has "transient" recovery, recovery time is defined as
        %   elapsed time between the end of the squeeze segment and the
        %   beginning of the first segment where the cell was recovered
        
        % "recovered" is when dI_rec comes within 8% of dI_np or higher
        rec_tol = 0.08;
        
        if num_rec_segs ~= 3
            warning('recovery time & category are hard-coded for devices with 3 recovery segments!');
        end
        
        rec_cat_description = '0 = instant, 1-2 = transient, 3 = prolonged';
        
        % cell was already recovered by the first recovery segment
        if (dI_np-dI_rec(1))/dI_np < rec_tol
            T_rec = (pulse_series(rec_k_start,1) - pulse_series(sq_k,2)) /Fs*N;
            rec_cat = 0;

        % cell didn't recover until the second recovery segment
        elseif (dI_np-dI_rec(2))/dI_np < rec_tol
            T_rec = (pulse_series(rec_k_start+1,1) - pulse_series(sq_k,2)) /Fs*N;
            rec_cat = 1;

        % cell didn't recover until the third recovery segment
        elseif (dI_np-dI_rec(3))/dI_np < rec_tol
            T_rec = (pulse_series(rec_k_start+2,1) - pulse_series(sq_k,2)) /Fs*N;
            rec_cat = 2;

        % by the third recovery segment, cell still hadn't recovered
        else
            T_rec = Inf;
            rec_cat = 3;

        end

        %% perform mNPS-r recovery curve fitting
        if fitflag

            if num_rec_segs ~= 3
                error('mNPS-r recovery fitting is hard-coded for devices with exactly 3 recovery segments');
            end

            % populate vector rT with time-points of recovery pulses
            rT = [ 0, ...
                (pulse_series(rec_k_start,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+1,2) - pulse_series(sq_k,1)) ] ...
                ./ Fs .* N; % in ms

            % populate vector rdI with current drop-amplitudes of recovery
            % pulses. Anticipate that rdI should increase
            %
            % NOTE: cannot approximate size for ellipsoid particle
            rdI = dI_rec';

            % use MATLAB fit() to fit a linear polynomial to recovery data.
            % fo has fields pertaining to fit parameters:
            %   fo.p1 is the slope of the line, not sign-bounded
            %   fo.p2 is the y-intercept, not sign-bounded
            % gof has fields pertaining to goodness of fit, but only
            % gof.rsquare will be used
            [fo, gof] = fit(rT',rdI','poly1');

        else
            fo.p1 = 0;
            fo.p2 = 0;
            gof.rsquare = 0;
        end

        out(k,:) = [start_index, I_baseline, dI_np, dI_c, dI_c_std, dT_np, dT_c, T_rec, fo.p1, fo.p2, gof.rsquare, rec_cat];
        out_cols = {'start_ix', 'I_baseline', 'dI_np', 'dI_c', 'dI_c_std', ...
            'dT_np', 'dT_c', 'T_rec', 'fo_p1', 'fo_p2', 'gof_rsquare', 'rec_cat'};
        out_units = {'index', 'data units', 'data units', 'data units', 'data units', ...
            'ms',    'ms',   'ms',    '',      '',      '',        'categorical'};

    end

    if fitflag

        Pix_SS = get(0,'screensize');

        figf = figure(43);
        title('Curve Fit'),
        ffunc = @(x) fo.p2+fo.p1*x;
        fplot(ffunc,[rT(1), rT(end)]),
        hold on, scatter(rT,rdI), hold off,
        figsize = [0.70 0.6 1/3 2/9]*Pix_SS(4);
        set(figf,'units','pixels','pos',figsize);

    end

    %% Section 9b: compute derived values

    calculated = zeros(size(out,1),9);

    % diameter (um)
    calculated(:,1) = ((out(:,3)./out(:,2)*De_np^2*L)./ ...
        (1+0.8*L/De_np*out(:,3)./out(:,2))).^(1/3);

    % strain (dimensionless)
    calculated(:,2) = (calculated(:,1)-wC)./calculated(:,1);

    % np velocity (mm/s = µm/ms)
    calculated(:,3) = mean( npL_ref ./ (dT_np_segs) );

    % sq velocity (mm/s = µm/ms)
    calculated(:,4) = sqL./out(:,7);

    % deformed diameter (um)
    calculated(:,5) = 0.01*(pi/4*wC)*(((out(:,4)./out(:,2)*De_c^2*L)./ ...
        (1+0.8*L/De_c*out(:,4)./out(:,2))).^(1/3)).^2;

    % wCDI (dimensionless)
    calculated(:,6) = calculated(:,4)./calculated(:,3) .* calculated(:,1)./ch_height;

    % recovery time (ms)
    calculated(:,7) = out(:,8);

    % recovery rate
    calculated(:,8) = fo.p1;

    % recovery category
    calculated(:,9) = out(:,12);
    
    % column names & units
    calc_cols = {'diameter', 'strain', 'V_np',          'V_c', 'def_diameter', ...
        'wCDI',         'rec_time', 'rec_rate', 'rec_cat'};
    calc_units = {'µm', 'dimensionless', 'mm/s = µm/ms', 'mm/s = µm/ms', 'µm', ...
        'dimensionless', 'ms',      '',         'categorical'};

    %% concatenate pulse data with calculated values
    
    OUT_array = [out(:,1:11), calculated]; % don't include rec_cat twice
    
    % set names & units for the columns
    column_names = [out_cols(1:11), calc_cols];
    column_units = [out_units(1:11), calc_units];
    
end

