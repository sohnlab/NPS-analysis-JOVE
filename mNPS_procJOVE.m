function output_table = mNPS_procJOVE(filepath, ch_height, De_np, wC, thresholds, sampleRate, ASLS_param, eventlength_filt)
% [ output_matrix ] = sNPS( start, number_of_files, thresholds )
%   Reads all sNPS data, analyzes data and returns final output matrix.
%   Needs a vector of 2 thresholds for initial thresholding. Afterwards,
%       QC, thresholding, and analysis for individual pulses will need to
%       proceed with user input.
% INPUTS:
%   filepath = char
%       the path to a .mat file containing the variable `data` (1xn double)
%   ch_height = double [um]
%       channel height measured from the SU-8 wafer
%   De_np = double [um]
%       effective diameter for the node-pore segments
%       measured by running calibration particles through the REF device
%   wC = double [um]
%       contraction channel width
%   thresholds (optional) = 1x2 vector of doubles
%       [low_threshold, high_threshold]
%       default = [1e-4, 1e-3]
%       pass [] or '' to use default thresholds
%   sampleRate (optional) = int [Hz]
%       default = 50e3
%   ASLS_param (optional) = struct
%       baseline fitting parameters passed to ASLS.m
%       if not provided, fills with default parameter values
%   eventlength_filt (optional)
%       expected # of filtered samples for a whole cell transit (default 2000)

    %% parse inputs

    if nargin < 5 || isempty(thresholds)
        thresholds = [1e-4, 1e-3];
        fprintf('Auto thresholds set to %3.2e, %3.2e\n',thresholds);
    end

    if nargin < 6 || isempty(sampleRate)
        sampleRate = 50000;
        fprintf('default sample rate used: %d Hz\n', sampleRate);
    end

    % default ASLS parameters
    if nargin<7 || isempty(ASLS_param)
        ASLS_param = struct();
        ASLS_param.lambda = 1e9; % default 1e5; larger=smoother, smaller=wiggly-er (may not be unit-independent)
        ASLS_param.p = 3e-3; % default 0.01; 0>p>1 (as low as possible while still converging)
        ASLS_param.noise_margin = 1e-4; % default 0; allows baseline to sit within the baseline noise
        ASLS_param.max_iter = 20; % default 5; just make sure it converges
    end

    if nargin<8 || isempty(eventlength_filt)
        eventlength_filt = 2000;
    end

    %% read all, search for pulses, get column information

    % load data
    load(filepath,'data');
    if ~( size(data,1)==1 && size(data,2)>2 )
        error('the variable `data` in filepath must be a 1xn array of doubles');
    end

    [all_out, ~, ~, outcols, outunits, rec_des] = ...
        mNPS_readJOVE(data, sampleRate, ch_height, De_np, wC, thresholds, false, false, ASLS_param);

    %% remove duplicate files to read

    [uni_win, output_matrix] = mNPS_cleanKim(all_out);

    clear all_out
    %% analyze all time windows, prompt user for input to see if the pulse looks good

    good_index = 1; % index through output of good pulses
    i = 1;
    new_th = thresholds; % reset threshold values
    warning('off', 'curvefit:cfit:subsasgn:coeffsClearingConfBounds'); % annoying warning when fit fails to converge

    while (i < size(uni_win,1))

        % checking for working pulses
        searchflag = true; % true to search, false to stop
        retryflag = true; % true for first attempt, false if trying again with default thresholds

        while (searchflag)
            try
                fprintf('\nNow reading index: %d\n',uni_win(i));
                fprintf('Progress: %2.1f %%\n',i/length(uni_win)*100);

                % set window size parameters
                min_windowsize_filt = ceil(sampleRate/20); % for LPF padding, ensure window size cover >1sec
                startoffset_filt = 200; % # filtered samples included before the first detected peak
                windowsize_filt = eventlength_filt * 1.2; % # filtered samples w/ 20% buffer
                % set window indices
                startix_filt = uni_win(i)-startoffset_filt;
                endix_filt = startix_filt + max(min_windowsize_filt,windowsize_filt);
                startix = max(1, 20*startix_filt); % ensure valid index
                endix = min(length(data), 20*endix_filt); % ensure valid index
                % extract iterdata window & check length
                iterdata = data(startix:endix);
                if length(iterdata) < sampleRate+1
                    warning('iterdata has only %u samples, but sampleRate=%u; this can cause problems with pad_data in LPF', length(iterdata), sampleRate);
                end

                % measure a new pulse
                [iter_out, emptyflag] = mNPS_readJOVE(iterdata, sampleRate, ch_height, De_np, wC, new_th, false, false, ASLS_param);
                % iter_out: output of one iteration
                % emptyflag: skip pulse if TRUE
                % auto: values for computing auto-threshold value
                % 
                % function arguments: read iterdata from uni_win entries
                % plot and fit off for speed


                if emptyflag
                    fprintf('No pulse! Skipping...\n');
                    i = i+1;
                    searchflag = true;
                else
                    [iter_out, ~, auto] = mNPS_readJOVE(iterdata, sampleRate, ch_height, De_np, wC, new_th, true, true, ASLS_param);
                    searchflag = false;
                end

            catch ME

                %% print errors to command line for debug purposes
                fprintf('-----\n%s\n',ME.identifier);
                for errorstack_i = 1:length(ME.stack)
                    fprintf('Line: %d --- %s\n',ME.stack(errorstack_i).line,ME.stack(errorstack_i).name);
                end
                if retryflag % retry once with default thresholds
                    new_th = thresholds;
                    searchflag = true;
                    fprintf('Error occured, retrying with default thresholds...\n');
                    retryflag = false;
                else % fail on second try -> skip entirely
                    new_th = thresholds;
                    i = i + 1; % skip this pulse
                    fprintf('Error occured, skipping this file...\n');
                    retryflag = true;
                end
            end

            if i > size(uni_win,1) % return if finished
                break
            end
        end

        % prompt for input to determine next operation
        fprintf(['Is anything wrong?\n' ...
                    'Press RETURN to skip this pulse\n' ...
                    'Enter X or any unrecognized character to stop analyzing\n\n' ...
                    'Enter P or . or / to save this data (auto-choose window)\n' ...
                    'Enter // to choose the window to save\n' ...
                    'Enter T to adjust the top threshold only\n' ...
                    'Enter Y to adjust the bottom threshold only\n' ...
                    'Enter + or = or '' to auto-adjust according to top threshold\n' ...
                    'Or, just enter a new top threshold to automatically set both\n' ...
                    ]);
        fprintf('\n%d cell(s) saved\n',good_index-1);
        fprintf('\nCurrent thresholds: %3.2e, %3.2e\n',new_th);
        OK = input('---\n','s');

        %% input case structures
        switch OK
            
            % empty input: skip the pulse
            case []

            fprintf('Skipping this pulse...\n');
            i = i+1;
            new_th = thresholds;

            % save data, auto-choose window
            case {'P','p', '.', '/'}

%             fprintf('Ok, displaying windows and file numbers...\n');

            % make and display table of pulses and indices
            indices = iter_out(iter_out(:,1) > 100);
            winds = 1:length(indices);
            table_data = [winds', indices];

            % clean up empty table entries
            cci = 1;
            stopc = size(table_data,1);
            while(cci <= stopc)
                if table_data(cci,2) == 0
                    table_data(cci,:) = [];
                    stopc = stopc - 1;
                else
                    cci = cci + 1;
                end
            end


            if ~isempty(table_data) % make sure WindowTable is NOT empty
    %             iter_out_index = input('Select window to save:\n---\n');

                disp( array2table(table_data,'VariableNames',{'Window','StartIndex'}) );
                iter_out_index = 1;
                fprintf('Ok, saving data...\n');

                if (0 < iter_out_index) && (iter_out_index <= max(winds))
                    output_matrix(good_index,:) = iter_out(iter_out_index,:); % save to output matrix
                    output_matrix(good_index,1) = output_matrix(good_index,1) + uni_win(i) - 200;
                    good_index = good_index + 1;
                    i = i+1;
                    new_th = thresholds; % reset thresholds
                else
                    fprintf('Unrecognized input.\n');
                    beep
                end

            else
                fprintf('Pulse table is empty; force retry\n');
                new_th = thresholds; % reset thresholds
            end

            % save data, user picks window
            case {'//'}

%             fprintf('Ok, displaying windows and file numbers...\n');

            % make and display table of pulses and indices
            indices = iter_out(:,1);
            winds = 1:length(indices);
            table_data = [winds', indices];

            % clean up empty table entries
            cci = 1;
            stopc = size(table_data,1);
            while(cci <= stopc)
                if table_data(cci,2) == 0
                    table_data(cci,:) = [];
                    stopc = stopc - 1;
                else
                    cci = cci + 1;
                end
            end

            if ~isempty(table_data) % make sure WindowTable is NOT empty

                disp( array2table(table_data,'VariableNames',{'Window','StartIndex'}) );
                iter_out_index = input('Select window to save:\n---\n');
                fprintf('Ok, saving data...\n');

                if (0 < iter_out_index) && (iter_out_index <= max(winds))
                    output_matrix(good_index,:) = iter_out(iter_out_index,:); % save to output matrix
                    output_matrix(good_index,1) = output_matrix(good_index,1) + uni_win(i) - 200;
                    good_index = good_index + 1;
                    i = i+1;
                    new_th = thresholds; % reset thresholds
                else
                    fprintf('Unrecognized input.\n');
                    beep
                end

            else
                fprintf('Pulse table is empty; force retry\n');
                new_th = thresholds; % reset thresholds
            end

            % in case of bad threshold, ask for new thresholds and retry

            % manually set new top threshold
            case {'T','t'}

            th_input = input('Please input new top threshold\n---\n');

                if isempty(th_input)
                    fprintf('Unrecognized input: stopping and dumping data.\n');
                    break

                else
                    new_th(2) = th_input;
                    fprintf('New top threshold: %3.2e\n',new_th(2));
                end

            % manually set new bottom threshold
            case {'Y','y'}

            th_input = input('Please input new bottom threshold\n---\n');

                if isempty(th_input)
                    fprintf('Unrecognized input: stopping and dumping data.\n');
                    break

                else
                    new_th(1) = th_input;
                    fprintf('New bottom threshold: %3.2e\n',new_th(1));
                end

            % automatically compute threshold values based on max vals from difference vector
            case {'+', '=', ''''}
                new_th(2) = auto*0.85;
                new_th(1) = new_th(2)*0.12;
                if new_th(1) < thresholds(1)
                    new_th(1) = thresholds(1);
                end
                fprintf('--------------\nSet new bottom threshold to %3.2e\n',new_th(1));
                fprintf('\nSet new upper threshold to %3.2e\n',new_th(2));

            otherwise

                % manually set new top threshold and automatically compute bottom threshold
                if length(OK) > 1
                    th_input = str2double(OK);
                    new_th(1) = th_input/10;
                    new_th(2) = th_input;
                    fprintf('New bottom threshold: %3.2e\n',new_th(1));
                    fprintf('New top threshold: %3.2e\n',new_th(2));

                % unrecognized character: stop analyzing data and save processed events so far
                else
                    fprintf('Unrecognized input: stopping and dumping data.\n');
                    break
                end
        end
    end

    %% finish up

    % remove empty rows from the output arrray
    output_data = output_matrix( ~all(output_matrix==0,2), :);
    
    % don't include recovery time twice
    output_data = output_data(:, ~strcmp(outcols, 'T_rec'));
    column_names = outcols(~strcmp(outcols, 'T_rec'));
    column_units = outunits(~strcmp(outcols, 'T_rec'));

    % description for the recovery category
    column_descriptions = cell(size(column_names));
    column_descriptions(:) = {''};
    column_descriptions(strcmp(column_names, 'rec_cat')) = {rec_des};
    
    % convert output array to a table with variable names & info
    output_table = array2table(output_data, 'variablenames', column_names);
    output_table.Properties.DimensionNames{2} = 'cell_data';
    output_table.Properties.VariableUnits = column_units;
    output_table.Properties.VariableDescriptions = column_descriptions;

    % play an alert sound to signal that processing is finished
    t = linspace(0,1,2^16); % time-samples
    Fs = 2^16; % sampling frequency
    y = 0.3*exp(-4*t).*sin(t*2*pi*440); % sinusoid
    sound(y,Fs);
    
    % remove duplicate cell detections (keep last-recorded instance)
    output_table = remove_duplicate_rows(output_table);

    fprintf('Done reading, check output!\n');

end

%% remove obvious duplicates (keep last-recorded instance)

function new_table = remove_duplicate_rows(old_table)
    [~,ia,~] = unique(old_table.start_ix, 'last');
    new_table = old_table(ia,:);
end


