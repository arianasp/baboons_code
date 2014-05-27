 function [ interactions ] = dyadic_interactions( xs, ys, day_start_idxs, day_range, a, b, noise_thresh )
%Finds and categorizes all of the dyadic interactions between two
%individuals, a and b, and outputs information about them
%INPUTS:
%   xs: [NxT matrix] of x positions of each individual over time
%   ys: [NxT matrix] of y positions of each individual over time
%   day_start_idxs: [n_days vector] of indices at which each day started
%   day_range: [variable length vector] of which days to include in the analysis
%   a: [number] index of the first individual
%   b: [number] index of the second individual
%   noise_thresh: [number] how much noise in inter-dyadic distance to allow before it is considered that the dyadic distance has left the current local minimum or maximum
%OUTPUTS:
%   interactions: [struct] containing information about all interactions. The struct contains the following fields:]
%       .day: day on which the interaction occurred
%       .noise_thresh: noise threshold used (how much the dyadic distance
%           had to change in order to leave the current local min or max)
%       .time_idxs: [t1 t2 t3] (times involved in the interaction)
%       .strength_add: additive "strength" of the interaction
%       .strength_mult: multiplicative "strength" of the interaction
%       .disp_add: additive "disparity" of the interaction
%       .disp_mult: multiplicative "displarity" of the interaction
%       .dyadic_dist_diffs: difference between dyadic distance between time
%           t1 and t2 and between t2 and t3
%       .travel_diffs: difference between the amount a and b traveled 
%           during the time between t1 and t2 and between t2 and t3 (total
%           displacement, not path length)
%       .type: type of interaction ('push','pull','repel','anchor'
%       .leader: which individual was the "leader" in the interaction (i.e.
%           the one who pushed, pulled, anchored, or repelled the other
%           individual)
%       .follower: which individual was the "follower" in the interaction
%           (i.e. the one who was pushed, pulled, anchored, or repelled by
%           the other individual)

T = size(xs,2);
day_start_idxs = [day_start_idxs T+1];

%get dyadic distance over time
[ dyadic_dists ] = dyadic_dist_over_time( xs, ys, a, b);


%struct to hold all interactions
interactions = {};

idx = 1;
for d = day_range
    
    t_range = day_start_idxs(d):(day_start_idxs(d+1)-1);
    dyadic_dists_curr = dyadic_dists(t_range);
    
    %break into sequences that contain no NaNs
    [seqs seq_starts] = get_non_nan_sequences(dyadic_dists_curr',3);
    
    %gets mins and maxes for each sequence and aggregate them together into
    %one list
    mins_maxes = [];
    for i = 1:length(seqs)
        mins_maxes_curr = find_max_min_points(seqs{i},noise_thresh);
        mins_maxes = [mins_maxes; mins_maxes_curr + seq_starts(i)-1];
    end
    
    %adjust the indices on mins_maxes
    mins_maxes = mins_maxes + day_start_idxs(d) - 1;
    
    %extract triads of times
    min_max_min = nan(length(mins_maxes),3);
    max_min_max = nan(length(mins_maxes),3);
    idx2 = 1;
    for i = 1:(size(mins_maxes,1)-1)
        min_max_min(idx2,1) = mins_maxes(i,1);
        min_max_min(idx2,2) = mins_maxes(i,2);
        min_max_min(idx2,3) = mins_maxes(i+1,1);
        max_min_max(idx2,1) = mins_maxes(i,2);
        max_min_max(idx2,2) = mins_maxes(i+1,1);
        max_min_max(idx2,3) = mins_maxes(i+1,2);
        idx2 = idx2 + 1;
    end
    
    %remove rows with nans
    min_max_min = remove_rows_with_any_nans( min_max_min );
    max_min_max = remove_rows_with_any_nans( max_min_max );
    
    %change in location of each individual between consecutive time points
    nrow_min_max_min = size(min_max_min,1);
    nrow_max_min_max = size(max_min_max,1);
    
    min_max_min_dx_a = diff(reshape(xs(a,min_max_min),nrow_min_max_min,3),1,2);
    min_max_min_dy_a = diff(reshape(ys(a,min_max_min),nrow_min_max_min,3),1,2);
    max_min_max_dx_a = diff(reshape(xs(a,max_min_max),nrow_max_min_max,3),1,2);
    max_min_max_dy_a = diff(reshape(ys(a,max_min_max),nrow_max_min_max,3),1,2);
    min_max_min_dx_b = diff(reshape(xs(b,min_max_min),nrow_min_max_min,3),1,2);
    min_max_min_dy_b = diff(reshape(ys(b,min_max_min),nrow_min_max_min,3),1,2);
    max_min_max_dx_b = diff(reshape(xs(b,max_min_max),nrow_max_min_max,3),1,2);
    max_min_max_dy_b = diff(reshape(ys(b,max_min_max),nrow_max_min_max,3),1,2);
    
    %distance traveled by each individual during consecutive time points
    min_max_min_dist_a = sqrt(min_max_min_dx_a.^2 + min_max_min_dy_a.^2);
    min_max_min_dist_b = sqrt(min_max_min_dx_b.^2 + min_max_min_dy_b.^2);
    max_min_max_dist_a = sqrt(max_min_max_dx_a.^2 + max_min_max_dy_a.^2);
    max_min_max_dist_b = sqrt(max_min_max_dx_b.^2 + max_min_max_dy_b.^2);
    
    %dyadic distances
    min_max_min_dyadic_dists = reshape(dyadic_dists(min_max_min),nrow_min_max_min,3);
    max_min_max_dyadic_dists = reshape(dyadic_dists(max_min_max),nrow_max_min_max,3);
    
    %change in dyadic distances between consecutive time points
    min_max_min_dyadic_dist_diffs = diff(min_max_min_dyadic_dists,1,2);
    max_min_max_dyadic_dist_diffs = diff(max_min_max_dyadic_dists,1,2);
    
    %vector between individual a and b during each time
    x_a_min_max_min = reshape(xs(a,min_max_min),nrow_min_max_min,3);
    y_a_min_max_min = reshape(ys(a,min_max_min),nrow_min_max_min,3);
    x_b_min_max_min = reshape(xs(b,min_max_min),nrow_min_max_min,3);
    y_b_min_max_min = reshape(ys(b,min_max_min),nrow_min_max_min,3);
    
    x_a_max_min_max = reshape(xs(a,max_min_max),nrow_max_min_max,3);
    y_a_max_min_max = reshape(ys(a,max_min_max),nrow_max_min_max,3);
    x_b_max_min_max = reshape(xs(b,max_min_max),nrow_max_min_max,3);
    y_b_max_min_max = reshape(ys(b,max_min_max),nrow_max_min_max,3);
    
    dx_ab_min_max_min = x_b_min_max_min - x_a_min_max_min;
    dy_ab_min_max_min = y_b_min_max_min - y_a_min_max_min;
    
    dx_ab_max_min_max = x_b_max_min_max - x_a_max_min_max;
    dy_ab_max_min_max = y_b_max_min_max - y_a_max_min_max;
    
    %get difference in distance traveled by a and during t1-->t2 (col 1)
    %and t2--> t3 (col 2) during each event
    min_max_min_travel_diffs = min_max_min_dist_a - min_max_min_dist_b;
    max_min_max_travel_diffs = max_min_max_dist_a - max_min_max_dist_b;
    %get sum of distance traveled by both individuals during t1-->t2 (col 1)
    %and t2-->t3 (col 2) during each event
    min_max_min_travel_sums = min_max_min_dist_a + min_max_min_dist_b;
    max_min_max_travel_sums = max_min_max_dist_a + max_min_max_dist_b;
    
    %get strengths
    min_max_min_strengths_add = (abs(min_max_min_dyadic_dist_diffs(:,1)) + abs(min_max_min_dyadic_dist_diffs(:,2)))./...
        ((min_max_min_dyadic_dists(:,1)+min_max_min_dyadic_dists(:,2))+(min_max_min_dyadic_dists(:,2)+min_max_min_dyadic_dists(:,3)));
    max_min_max_strengths_add = (abs(max_min_max_dyadic_dist_diffs(:,1)) + abs(max_min_max_dyadic_dist_diffs(:,2)))./...
        ((max_min_max_dyadic_dists(:,1)+max_min_max_dyadic_dists(:,2))+(max_min_max_dyadic_dists(:,2)+max_min_max_dyadic_dists(:,3)));
    min_max_min_strengths_mult = (abs(min_max_min_dyadic_dist_diffs(:,1)) .* abs(min_max_min_dyadic_dist_diffs(:,2)))./...
        ((min_max_min_dyadic_dists(:,1)+min_max_min_dyadic_dists(:,2)).*(min_max_min_dyadic_dists(:,2)+min_max_min_dyadic_dists(:,3)));
    max_min_max_strengths_mult = (abs(max_min_max_dyadic_dist_diffs(:,1)) .* abs(max_min_max_dyadic_dist_diffs(:,2)))./...
        ((max_min_max_dyadic_dists(:,1)+max_min_max_dyadic_dists(:,2)).*(max_min_max_dyadic_dists(:,2)+max_min_max_dyadic_dists(:,3)));
    
    %get disparities
    min_max_min_disps_add = (abs(min_max_min_travel_diffs(:,1)) + abs(min_max_min_travel_diffs(:,2)))./...
        (min_max_min_travel_sums(:,1) + min_max_min_travel_sums(:,2));
    max_min_max_disps_add = (abs(max_min_max_travel_diffs(:,1)) + abs(max_min_max_travel_diffs(:,2)))./...
        (max_min_max_travel_sums(:,1) + max_min_max_travel_sums(:,2));
    min_max_min_disps_mult = (abs(min_max_min_travel_diffs(:,1)) .* abs(min_max_min_travel_diffs(:,2)))./...
        (min_max_min_travel_sums(:,1) .* min_max_min_travel_sums(:,2));
    max_min_max_disps_mult = (abs(max_min_max_travel_diffs(:,1)) .* abs(max_min_max_travel_diffs(:,2)))./...
        (max_min_max_travel_sums(:,1) .* max_min_max_travel_sums(:,2));
    
    for i = 1:nrow_min_max_min
        interactions(idx).day = d;
        interactions(idx).noise_thresh = noise_thresh;
        interactions(idx).time_idxs = min_max_min(i,:);
        interactions(idx).strength_add = min_max_min_strengths_add(i);
        interactions(idx).strength_mult = min_max_min_strengths_mult(i);
        interactions(idx).disp_add = min_max_min_disps_add(i);
        interactions(idx).disp_mult = min_max_min_disps_mult(i);
        interactions(idx).dyadic_dist_diffs = min_max_min_dyadic_dists(i,:);
        interactions(idx).travel_diffs = min_max_min_travel_diffs(i,:);
        if min_max_min_travel_diffs(i,1) > 0
            if min_max_min_travel_diffs(i,2) > 0
                interactions(idx).type = 'anchor';
                interactions(idx).leader = b;
                interactions(idx).follower = a;
                interactions(idx).leader_disp_12 = [min_max_min_dx_b(i,1); min_max_min_dy_b(i,1)];
                interactions(idx).leader_disp_23 = [min_max_min_dx_b(i,2); min_max_min_dy_b(i,2)];
                interactions(idx).follower_disp_12 = [min_max_min_dx_a(i,1); min_max_min_dy_a(i,1)];
                interactions(idx).follower_disp_23 = [min_max_min_dx_a(i,2); min_max_min_dy_a(i,2)];
                interactions(idx).dyad_vector_t2 = [dx_ab_min_max_min(i,2); dy_ab_min_max_min(i,2)];
            elseif min_max_min_travel_diffs(i,2) < 0
                interactions(idx).type = 'pull';
                interactions(idx).leader = a;
                interactions(idx).follower = b;
                interactions(idx).leader_disp_12 = [min_max_min_dx_a(i,1); min_max_min_dy_a(i,1)];
                interactions(idx).leader_disp_23 = [min_max_min_dx_a(i,2); min_max_min_dy_a(i,2)];
                interactions(idx).follower_disp_12 = [min_max_min_dx_b(i,1); min_max_min_dy_b(i,1)];
                interactions(idx).follower_disp_23 = [min_max_min_dx_b(i,2); min_max_min_dy_b(i,2)];
                interactions(idx).dyad_vector_t2 = [-dx_ab_min_max_min(i,2); -dy_ab_min_max_min(i,2)];
            end
        elseif min_max_min_travel_diffs(i,1) < 0
            if min_max_min_travel_diffs(i,2) > 0
                interactions(idx).type = 'pull';
                interactions(idx).leader = b;
                interactions(idx).follower = a;
                interactions(idx).leader_disp_12 = [min_max_min_dx_b(i,1); min_max_min_dy_b(i,1)];
                interactions(idx).leader_disp_23 = [min_max_min_dx_b(i,2); min_max_min_dy_b(i,2)];
                interactions(idx).follower_disp_12 = [min_max_min_dx_a(i,1); min_max_min_dy_a(i,1)];
                interactions(idx).follower_disp_23 = [min_max_min_dx_a(i,2); min_max_min_dy_a(i,2)];
                interactions(idx).dyad_vector_t2 = [dx_ab_min_max_min(i,2); dy_ab_min_max_min(i,2)];
            elseif min_max_min_travel_diffs(i,2) < 0
                interactions(idx).type = 'anchor';
                interactions(idx).leader = a;
                interactions(idx).follower = b;
                interactions(idx).leader_disp_12 = [min_max_min_dx_a(i,1); min_max_min_dy_a(i,1)];
                interactions(idx).leader_disp_23 = [min_max_min_dx_a(i,2); min_max_min_dy_a(i,2)];
                interactions(idx).follower_disp_12 = [min_max_min_dx_b(i,1); min_max_min_dy_b(i,1)];
                interactions(idx).follower_disp_23 = [min_max_min_dx_b(i,2); min_max_min_dy_b(i,2)];
                interactions(idx).dyad_vector_t2 = [-dx_ab_min_max_min(i,2); -dy_ab_min_max_min(i,2)];
            end
        end

        idx = idx + 1;
        
    end
    
    for i = 1:nrow_max_min_max
        interactions(idx).day = d;
        interactions(idx).noise_thresh = noise_thresh;
        interactions(idx).time_idxs = max_min_max(i,:);
        interactions(idx).strength_add = max_min_max_strengths_add(i);
        interactions(idx).strength_mult = max_min_max_strengths_mult(i);
        interactions(idx).disp_add = max_min_max_disps_add(i);
        interactions(idx).disp_mult = max_min_max_disps_mult(i);
        interactions(idx).dyadic_dist_diffs = max_min_max_dyadic_dists(i,:);
        interactions(idx).travel_diffs = max_min_max_travel_diffs(i,:);
        if max_min_max_travel_diffs(i,1) > 0
            if max_min_max_travel_diffs(i,2) > 0
                interactions(idx).type = 'repel';
                interactions(idx).leader = b;
                interactions(idx).follower = a;
                interactions(idx).leader_disp_12 = [max_min_max_dx_b(i,1); max_min_max_dy_b(i,1)];
                interactions(idx).leader_disp_23 = [max_min_max_dx_b(i,2); max_min_max_dy_b(i,2)];
                interactions(idx).follower_disp_12 = [max_min_max_dx_a(i,1); max_min_max_dy_a(i,1)];
                interactions(idx).follower_disp_23 = [max_min_max_dx_a(i,2); max_min_max_dy_a(i,2)];
                interactions(idx).dyad_vector_t2 = [dx_ab_max_min_max(i,2); dy_ab_max_min_max(i,2)];
            elseif max_min_max_travel_diffs(i,2) < 0
                interactions(idx).type = 'push';
                interactions(idx).leader = a;
                interactions(idx).follower = b;
                interactions(idx).leader_disp_12 = [max_min_max_dx_a(i,1); max_min_max_dy_a(i,1)];
                interactions(idx).leader_disp_23 = [max_min_max_dx_a(i,2); max_min_max_dy_a(i,2)];
                interactions(idx).follower_disp_12 = [max_min_max_dx_b(i,1); max_min_max_dy_b(i,1)];
                interactions(idx).follower_disp_23 = [max_min_max_dx_b(i,2); max_min_max_dy_b(i,2)];
                interactions(idx).dyad_vector_t2 = [-dx_ab_max_min_max(i,2); -dy_ab_max_min_max(i,2)];
            end
        elseif max_min_max_travel_diffs(i,1) < 0
            if max_min_max_travel_diffs(i,2) > 0
                interactions(idx).type = 'push';
                interactions(idx).leader = b;
                interactions(idx).follower = a;
                interactions(idx).leader_disp_12 = [max_min_max_dx_b(i,1); max_min_max_dy_b(i,1)];
                interactions(idx).leader_disp_23 = [max_min_max_dx_b(i,2); max_min_max_dy_b(i,2)];
                interactions(idx).follower_disp_12 = [max_min_max_dx_a(i,1); max_min_max_dy_a(i,1)];
                interactions(idx).follower_disp_23 = [max_min_max_dx_a(i,2); max_min_max_dy_a(i,2)];
                interactions(idx).dyad_vector_t2 = [dx_ab_max_min_max(i,2); dy_ab_max_min_max(i,2)];
            elseif max_min_max_travel_diffs(i,2) < 0
                interactions(idx).type = 'repel';
                interactions(idx).leader = a;
                interactions(idx).follower = b;
                interactions(idx).leader_disp_12 = [max_min_max_dx_a(i,1); max_min_max_dy_a(i,1)];
                interactions(idx).leader_disp_23 = [max_min_max_dx_a(i,2); max_min_max_dy_a(i,2)];
                interactions(idx).follower_disp_12 = [max_min_max_dx_b(i,1); max_min_max_dy_b(i,1)];
                interactions(idx).follower_disp_23 = [max_min_max_dx_b(i,2); max_min_max_dy_b(i,2)];
                interactions(idx).dyad_vector_t2 = [-dx_ab_max_min_max(i,2); -dy_ab_max_min_max(i,2)];
            end
        end
        
        idx = idx + 1;
        
    end








end

