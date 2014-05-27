%compute dyadic association for all pairs of individuals

%% Parameters
nbins = [20 20];
mins = [-100 -100];
maxes = [100 100];
dir = '/home/arianasp/baboons';
ndays = 14;
min_dist_sleep_site = 100;
n_randomizations = 100;
min_times = 500000; %minimum number of tracked seconds an individual must have to be included in the analysis

paras.nbins = nbins;
paras.mins = mins;
paras.maxes = maxes;
paras.ndays = ndays;
paras.n_randomizations = n_randomizations;
paras.min_dist_sleep_site = min_dist_sleep_site;
paras.min_times = min_times;

%% Load data
addpath(genpath(dir))
disp('loading data...')
load([dir '/data/matlab_processed/individual_level_metrics/pos_rel_to_centroid_level1.mat'],'xs_rel','ys_rel')
load([dir '/data/matlab_raw/day_start_idxs.mat'])
load([dir '/data/matlab_processed/individual_level_metrics/sleep_sites_level1.mat'],'troop_dist_from_sleep_site')
load([dir '/data/matlab_raw/baboon_info.mat'],'baboon_info')

xs = xs_rel(:,1:(day_start_idxs(ndays+1)-1));
ys = ys_rel(:,1:(day_start_idxs(ndays+1)-1));

%% Get only the individuals with enough data
disp('pre-processing data...')
times_tracked = sum(not(isnan(xs)),2);
inds_to_use = find(times_tracked >= paras.min_times);
xs = xs(inds_to_use,:);
ys = ys(inds_to_use,:);
baboon_info = baboon_info(inds_to_use);
paras.inds_used = inds_to_use;
N = length(inds_to_use);
T = size(xs,2);

%% Remove data from near the sleep site
troop_dist_from_sleep_site = troop_dist_from_sleep_site(1:(day_start_idxs(ndays+1)-1));
times_to_remove = find(troop_dist_from_sleep_site < min_dist_sleep_site);
xs(:,times_to_remove) = NaN;
ys(:,times_to_remove) = NaN;

%% Compute dyadic association for all pairs
disp('computing dyadic associations...')
assoc_data = cell(N,N);
assoc_mat = nan(N,N);
for i = 1:N
    i
    for j = i:N
        X = [xs(i,:)' ys(i,:)'];
        Y = [xs(j,:)' ys(j,:)'];
        dataij = dyadic_association(X,Y,mins,maxes,nbins);
        assoc_data{i,j} = dataij;
        assoc_data{j,i} = dataij;
        assoc_mat(i,j) = dataij.A_tot;
        assoc_mat(j,i) = dataij.A_tot;
    end
end

%% Randomize individual identities by day and recomputer dyadic associations
disp('computing randomized dyadic associations...')
assoc_data_rand = cell(N,N,n_randomizations);
assoc_mat_rand = nan(N,N,n_randomizations);
xy = nan(N,T,2);
xy(:,:,1) = xs;
xy(:,:,2) = ys;
for r = 1:n_randomizations
    [ xy_shuff ] = shuffle_ids_by_day( xy, day_start_idxs(1:ndays) );
    for i = 1:N
        disp([r i])
        for j = i:N
            X = [xy_shuff(i,:,1)' xy_shuff(i,:,2)'];
            Y = [xy_shuff(j,:,1)' xy_shuff(j,:,2)'];
            dataij = dyadic_association(X,Y,mins,maxes,nbins);
            assoc_data_rand{i,j,r} = dataij;
            assoc_data_rand{j,i,r} = dataij;
            assoc_mat_rand(i,j,r) = dataij.A_tot;
            assoc_mat_rand(j,i,r) = dataij.A_tot;
        end
    end
end

disp('saving data...')
save([dir '/positioning_analysis/dyadic_association/dyad_assoc_days1-' num2str(ndays) '_20_bins.mat'],'assoc_data','assoc_mat','assoc_data_rand','assoc_mat_rand','paras')
        
disp('done!')