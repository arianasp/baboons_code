%get all approach avoid interactions for each pair and each day
%put them all into a matrix called AA
%AA(i,j,d) = # of times i approached and j avoided on day d

%load data
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat');

xs_to_use = xs; %which x data to use
ys_to_use = ys; %which y data to use

n_inds = size(xs_to_use,1); % number of individuals
n_days = 14; %number of days in the data set

%create parameter struct
paras = {};
paras.prior_dist = 3;
paras.approach_dist = 2;
paras.final_dist = 3;
paras.approach_mvmt = 3;
paras.max_prior_avoider_mvmt = 1.5;
paras.max_subseq_approacher_mvmt = 1.5;
paras.min_avoid_mvmt = 3;
paras.min_time_between_interactions = 10;

AA = zeros(n_inds,n_inds,n_days); %matrix to hold approach-avoid interactions on each day

for d = 1:n_days %for each day
    
    d
    
    %get time range for that day
    t_start = day_start_idxs(d);
    t_end = day_start_idxs(d+1)-1;
    time_range = [t_start t_end];
    
    %for each pair of individuals
    for i = 1:n_inds
        for j = 1:n_inds
            if i ~= j %if they are the same individual, don't compute
                
                [ n_approach_avoids ] = get_approach_avoid_interactions( xs_to_use, ys_to_use, i, j, time_range, paras);
                AA(i,j,d) = n_approach_avoids;
            end
        end
    end
end

data = {};
data.paras = paras;
data.approach_avoids_by_day = AA;
data.approach_avoids_total = nansum(AA,3);

save('/Users/arianasp/Desktop/Dropbox/baboons_shared/network_data/dominance/dominance_network.mat','data')