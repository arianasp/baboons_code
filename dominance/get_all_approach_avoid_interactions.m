%get all approach avoid interactions for each pair and each day
%put them all into a matrix called AA
%AA(i,j,d) = # of times i approached and j avoided on day d

%load data
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat');

n_inds = size(xs,1); % number of individuals
n_days = 14; %number of days in the data set

%create a data struct to hold all data
data = {};

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
paras.days = 1:n_days;

data.paras = paras;

times = [];
approachers = [];
avoiders = [];

dom_mat = zeros(n_inds,n_inds);

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
                
                [ times_ij ] = get_approach_avoid_interactions( xs, ys, i, j, time_range, paras);
                times = [times times_ij];
                approachers = [approachers ones(1,length(times_ij))*i];
                avoiders = [avoiders ones(1,length(times_ij))*j];
                dom_mat(i,j) = dom_mat(i,j) + length(times_ij);
            end
        end
    end
end

data.dom_mat = dom_mat;
data.times = times;
data.approachers = approachers;
data.avoiders = avoiders;

