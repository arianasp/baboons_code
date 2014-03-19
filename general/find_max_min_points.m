function [ mins_maxes ] = find_max_min_points( data, noise_thresh )
%Finds the local maxes and mins of the 1D vector data, not counting
%fluctuations within a certain noise threshold (noise_thresh). This
%function is called in the function dyadic_interactions.m
%INPUTS:
%   data: [Tx1] vector of time series values
%   noise_thresh: threshold you have to surpass in order to be considered a
%           change in whether the data are going up or down
%OUTPUTS:
%   mins_maxes: [n_peaks x 2] vector giving the indexes to consecutive min 
%           (col 1) and max(col 2) values
%Note about NaNs: the function removes leading and trailing nans, but 
%cannot handle nans in the middle of data sequences (it will return an 
%empty matrix)

%remove leading and trailing nans
non_nans = find(not(isnan(data)));
if isempty(non_nans)
    mins_maxes = [];
    return
end
first = non_nans(1);
last = non_nans(end);

data = data(first:last);

if any(isnan(data))
    mins_maxes = [];
    return
    error('remove all but leading and trailing NaNs from data prior to entering it into this function')
end

T = length(data);

mins_maxes = nan(T,2);

%find out whether initial point is a local max or min
found = 0;
init_val = data(1);
i = 2;
while not(found) && (i+1) < length(data)
    curr_val = data(i);
    if abs(curr_val - init_val) > noise_thresh
        if curr_val - init_val > 0
            direction = 1;
            found = 1;
            %find minimum in preceding interval
            [val start_idx] = min(data(1:i));
            mins_maxes(1,1) = start_idx + first - 1;
        elseif curr_val - init_val < 0 
            direction = 0;
            found = 1;
            %find maximum in preceding interval
            [val start_idx] = max(data(1:i));
            mins_maxes(1,2) = start_idx + first - 1;
        end
    else
        i = i +1;
    end
end

if not(found)
    mins_maxes = [];
    return
end

if direction == 1
    idx = 1;
else
    idx = 2;
end

ref_idx = start_idx;
ref_val = data(start_idx);
for t = start_idx:T
    curr_val = data(t);
    if direction == 1 %if going up
        if curr_val > ref_val %if the current value is greater than the reference value
            ref_val = curr_val; %change the reference value to the current (greater) value
            ref_idx = t; %change the reference index accordingly
        elseif (curr_val - ref_val) < -noise_thresh %if curr val is lower than the ref value (beyond a noise thresh)
            mins_maxes(idx,2) = ref_idx+first-1; %record the local maximum (the last reference idx)
            direction = 0; %switch direction (to going down)
            idx = idx + 1; %increment the index
            ref_idx = t; %change the reference index to the previous value
            ref_val = data(ref_idx); %chance the reference value
        end
    elseif direction == 0 %if going down
       if curr_val < ref_val %if the current value is less than the reference value
           ref_val = curr_val; %change the reference value to the current (smaller) value
           ref_idx = t; %change the reference index accordingly
       elseif (curr_val - ref_val) > noise_thresh %if curr val is higher than the ref value (beyond a noise thresh)
           mins_maxes(idx,1) = ref_idx+first-1; %record the local minimum (the last reference idx)
           direction = 1; %switch direction (to going up)
           ref_idx = t; %change the reference index to the previous value
           ref_val = data(ref_idx); %change the reference value
       end
    end
            
end

%remove trailing NaNs
mins_maxes = mins_maxes(1:idx,:);
if isnan(mins_maxes(end,1))
    mins_maxes = mins_maxes(1:(end-1),:);
end
    


end

