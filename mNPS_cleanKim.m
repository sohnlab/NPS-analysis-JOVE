function [uni_win, output_matrix] = mNPS_cleanKim(all_out)
%mNPS_clean is used to remove errors from mNPS_read matrices

stop = size(all_out,1); % initialize stop to the number of pulses found
i = 1;

% erase stuff way too short
while (i < stop)
   if (all_out(i,7) < 2) % if it's faster than 2 ms
       all_out(i,:) = [];
       stop = stop - 1;      
   else
       i = i + 1;
   end
end

uni_win = all_out(:,1); % vector of unique samples to read

% clean up double counts
cci = 1;
stopc = size(uni_win,1);
while(cci < stopc)
    if (uni_win(cci+1) - uni_win(cci)) < 500
        uni_win(cci+1) = [];
        stopc = stopc - 1;
    else
        cci = cci + 1;
    end
end

output_matrix = zeros(length(uni_win),size(all_out,2)); % pre-allocate output matrix

end

