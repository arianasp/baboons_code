function [ times ] = get_approach_avoid_interactions( xs, ys, i, j, time_range, paras)
%Finds all the approach-avoid interactions between individual i and j
%during the time range given by time_range. i is the approacher, j is the
%avoider.
%NOTE: There are currently some arbitrary thresholds in this function. I
%just hard-coded them in for now because there were so many and I got tired
%of naming variables, but at some point I may want to generalize this
%function.
%INPUTS:
%   xs, ys: locations at all time poitns
%   i: the approacher
%   j: the avoider
%   time_range: a vector telling [first_time_idx last_time_idx]
%   paras: a struct containing various parameters, described below.
%An interaction counts as i approaching and j avoiding if (default values
%in parentheses):
%1. Distance between i and j > prior_dist (3) meters 10s prior to now
%2. i moved more than approach_mvmt (3) meters in past 10s
%3. i and j are less than approach_dist (2) meters apart now
%4. j moved less than max_prior_avoider_mvmt (1.5) meters in the past 30s
%5. i moved less than max_subseq_approacher_mvmt (1.5) meters in the next 30s
%6. j moved more than avoid_dist (3) meters in the next 10s
%7. Distance between i and j >  final_dist (3) meters 10s after now

%WARNING: DO NOT INSERT TIME RANGES THAT CROSS DAYS!

%Default values for parameters
if nargin == 5
    paras = {};
    paras.prior_dist = 3;
    paras.approach_dist = 2;
    paras.final_dist = 3;
    paras.approach_mvmt = 3;
    paras.max_prior_avoider_mvmt = 1.5;
    paras.max_subseq_approacher_mvmt = 1.5;
    paras.min_avoid_mvmt = 3;
    paras.min_time_between_interactions = 1;
end

%if time_range is in the wrong format, throw an error
if length(time_range) ~= 2
    error('time range is not in [start end] format')
end
    
n_approach_avoids = 0; %number of approach avoid interactions found
times = [];

idx = 1;
t = time_range(1)+30;
while t <= (time_range(2)-30)
    found = 0; %set 'found' flag to false
    %get current separation between i and j
    d_ij_curr = sqrt((xs(i,t)-xs(j,t))^2 + (ys(i,t)-ys(j,t))^2);
    if d_ij_curr < paras.approach_dist %if current distance apart is less than approach_dist
        d_ij_before = sqrt((xs(i,t-10)-xs(j,t-10))^2 + (ys(i,t-10) - ys(j,t-10))^2);
        if d_ij_before > paras.prior_dist %if distance 10 sec ago was greater than prior_dist
            d_ij_after = sqrt((xs(i,t+10)-xs(j,t+10))^2+(ys(i,t+10)-ys(j,t+10))^2);
            if d_ij_after > paras.final_dist %if distance apart 10 sec in future is greater than final_dist
                disp_i_past = nanmax(sqrt((xs(i,t)*ones(1,10)-xs(i,(t-10):(t-1))).^2 + (ys(i,t)*ones(1,10)-ys(i,(t-10):(t-1))).^2));
                if disp_i_past > paras.approach_mvmt %if i moved more than approach_mvmt in the past 10 sec
                    disp_j_past = nanmax(sqrt((xs(j,t)*ones(1,30)-xs(j,(t-30):(t-1))).^2+(ys(j,t)*ones(1,30)-ys(j,(t-30):(t-1))).^2));
                    if disp_j_past < paras.max_prior_avoider_mvmt % if j moved less than max_prior_avoider_mvmt in the past 30 sec
                        disp_i_fut = nanmax(sqrt((xs(i,(t+1):(t+30))-xs(i,t)*ones(1,30)).^2+(ys(i,(t+1):(t+30))-ys(i,t)*ones(1,30)).^2));
                        if disp_i_fut < paras.max_subseq_approacher_mvmt %if i moved less than max_subseq_approacher_mvmt in the next 30 sec
                            disp_j_fut = nanmax(sqrt((xs(j,(t+1):(t+10))-xs(j,t)*ones(1,10)).^2+(ys(j,(t+1):(t+10))-ys(j,t)*ones(1,10)).^2));
                            if disp_j_fut > paras.min_avoid_mvmt %if j moved more than min_avoid_mvmt in the next 10 sec
                                n_approach_avoids = n_approach_avoids + 1; %add an approach avoid interaction!
                                found = 1; %set found flag to true
                                times(idx) = t;
                            end
                        end
                    end
                end
            end
        end
    end
    
    %increment time. if found, increment by min_time_between_interactions
    if found 
        t = t + paras.min_time_between_interactions;
    else
        t = t + 1;
    end
    
end



end

