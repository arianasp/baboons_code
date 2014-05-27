function [ seqs, starts ] = get_non_nan_sequences( data, min_seq_length )
%Finds sequences of data within the array 'data' that do not have NaNs
%INPUTS:
%   data: [T x k array] of data, where there are T samples of
%   k-dimensional data
%   min_seq_length: mininum length of sequence to include
%OUTPUTS:
%   seqs: [cell array] with connected sequences that contain no NaNs
%   starts: [n_seqs x 1 vector] with the indexes to the start points of
%   each sequence

T = size(data,1);
k = size(data,2);

not_nans = sum(not(isnan(data)),2);
not_nans = not_nans == k;

seqs = {};
starts = [];
idx = 1;
start = 1;
for i = 1:T
    if i == T && i - start >= min_seq_length
        if not_nans(i)
            seqs{idx} = data(start:i,:);
            starts = [starts start];
        else
            seqs{idx} = data(start:(i-1),:);
            starts = [starts start];
        end
    elseif not_nans(i)
        continue
    else
        if i - start >= min_seq_length
            seqs{idx} = data(start:(i-1),:);
            starts = [starts start];
            idx = idx + 1;
        end
        start = i+1;
    end
end


end

