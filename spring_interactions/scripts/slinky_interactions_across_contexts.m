%slinky interaction networks across different contexts

%% Parameters

last_day = 14;
event_type = 'anchor';
min_dist_from_sleep_site = 100;
noise_thresh = 3;
min_strength = 0.1;
min_disp = 0.1;
figdir = '/Users/arianasp/Desktop/Baboons/output/push_pull/plots/interactions_across_contexts';
savedir = '/Users/arianasp/Desktop/Baboons/output/push_pull/networks_data';


min_speed_for_low_speed_travel = 0.1;
min_speed_for_high_speed_travel = 0.5;
min_eccentricity_in_progressions = 0.9;
max_tilt_in_progressions = 0.1;


%% Load data
load(['/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data/dyad_interactions_noise_thresh_' num2str(noise_thresh) '.mat'])
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/day_start_idxs.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/baboon_info.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/troop_level_metrics/order_params_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/terrain_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/speeds_step_10_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/sleep_sites_level1.mat')

N = size(interactions_all,1);
fig_subdir=[figdir '/' event_type '_noise_thresh_' num2str(noise_thresh) '_strength_' num2str(min_strength) '_disp_' num2str(min_disp)];
mkdir(fig_subdir) 

%% Get different time segments for different situations

%all times (within the specified day range and not including sleep sites) -
%this is the full data set to start with
all_times = day_start_idxs(1):(day_start_idxs(last_day+1)-1);
all_times = intersect(all_times, find(troop_dist_from_sleep_site > min_dist_from_sleep_site));

%very low speed (stationary) times
stationary_times = find(troop_speeds < min_speed_for_low_speed_travel);

%low speed travel times
low_speed_travel_times = find(troop_speeds >= min_speed_for_low_speed_travel .* troop_speeds < min_speed_for_high_speed_travel);

%high speed travel times
high_speed_travel_times = find(troop_speeds >= min_speed_for_high_speed_travel);

%progression times
progression_times = find(order_params.eccentricity > min_eccentricity_in_progressions .* order_params.tilt < max_tilt_in_progressions);
progression_times = intersect(progression_times, high_speed_travel_times);

%high speed non-progression times
non_progression_times = find(not(order_params.eccentricity > min_eccentricity_in_progressions .* order_params.tilt < max_tilt_in_progressions));
non_progression_times = intersect(non_progression_times, high_speed_travel_times);

%% get event matrices for each situation
disp('getting event matrices...')

[ events_all, direc_all ] = event_count_matrix( interactions_all, event_type, 'mult', min_strength, min_disp, all_times );
[ events_stationary, direc_stationary ] = event_count_matrix( interactions_all, event_type, 'mult', min_strength, min_disp, stationary_times );
[ events_low_speed_travel, direc_low_speed_travel ] = event_count_matrix( interactions_all, event_type, 'mult', min_strength, min_disp, low_speed_travel_times );
[ events_high_speed_travel, direc_high_speed_travel ] = event_count_matrix( interactions_all, event_type, 'mult', min_strength, min_disp, high_speed_travel_times );
[ events_progressions, direc_progressions ] = event_count_matrix( interactions_all, event_type, 'mult', min_strength, min_disp, progression_times );
[ events_non_progressions, direc_non_progressions ] = event_count_matrix( interactions_all, event_type, 'mult', min_strength, min_disp, non_progression_times );


%% remove individual 16 from everything
if length(baboon_info) == 26
    events_all = events_all([1:15 17:N],[1:15,17:N]); direc_all = direc_all([1:15 17:N],[1:15,17:N]);
    events_stationary = events_stationary([1:15 17:N],[1:15,17:N]); direc_stationary = direc_stationary([1:15 17:N],[1:15,17:N]);
    events_low_speed_travel = events_low_speed_travel([1:15 17:N],[1:15,17:N]); direc_low_speed_travel = direc_low_speed_travel([1:15 17:N],[1:15,17:N]);
    events_high_speed_travel = events_high_speed_travel([1:15 17:N],[1:15,17:N]); direc_high_speed_travel = direc_high_speed_travel([1:15 17:N],[1:15,17:N]);
    events_progressions = events_progressions([1:15 17:N],[1:15,17:N]); direc_progressions = direc_progressions([1:15 17:N],[1:15,17:N]);
    events_non_progressions = events_non_progressions([1:15 17:N],[1:15,17:N]); direc_non_progressions = direc_non_progressions([1:15 17:N],[1:15,17:N]);
    baboon_info = [baboon_info(1:15) baboon_info(17:26)];
    N = 25;
end


%% Rank matrices
disp('ranking event matrices...')
[ new_mat, ranks_all ] = rank_mat( direc_all, 'binary_mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end
new_mat(find(isnan(new_mat)))=0;
direc_all_ranked = new_mat;

[ new_mat, ranks_stationary ] = rank_mat( direc_stationary, 'binary_mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end
new_mat(find(isnan(new_mat)))=0;
direc_stationary_ranked = new_mat;

[ new_mat, ranks_low_speed_travel ] = rank_mat( direc_low_speed_travel, 'binary_mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end
new_mat(find(isnan(new_mat)))=0;
direc_low_speed_travel_ranked = new_mat;

[ new_mat, ranks_high_speed_travel ] = rank_mat( direc_high_speed_travel, 'binary_mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end
new_mat(find(isnan(new_mat)))=0;
direc_high_speed_travel_ranked = new_mat;

[ new_mat, ranks_progressions ] = rank_mat( direc_progressions, 'binary_mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end
new_mat(find(isnan(new_mat)))=0;
direc_progressions_ranked = new_mat;

[ new_mat, ranks_non_progressions ] = rank_mat( direc_non_progressions, 'binary_mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end
new_mat(find(isnan(new_mat)))=0;
direc_non_progressions_ranked = new_mat;


%% Store data in struct
network_data = {};
network_data.paras = {};
network_data.paras.event_type = event_type;
network_data.paras.day_range = 1:last_day;
network_data.paras.min_dist_from_sleep_site = min_dist_from_sleep_site;
network_data.paras.noise_thresh = noise_thresh;
network_data.paras.min_strength = min_strength;
network_data.paras.min_disp = min_disp;
network_data.paras.min_speed_for_low_speed_travel = min_speed_for_low_speed_travel;
network_data.paras.min_speed_for_high_speed_travel = min_speed_for_high_speed_travel;
network_data.paras.min_eccentricity_in_progressions = min_eccentricity_in_progressions;
network_data.paras.max_tilt_in_progressions = max_tilt_in_progressions;
network_data.baboon_info = baboon_info;
network_data.events_all = events_all;
network_data.events_stationary = events_stationary;
network_data.events_low_speed_travel = events_low_speed_travel;
network_data.events_high_speed_travel = events_high_speed_travel;
network_data.events_progressions = events_progressions;
network_data.events_non_progressions = events_non_progressions;
network_data.direc_all = direc_all;
network_data.direc_stationary = direc_stationary;
network_data.direc_low_speed_travel = direc_low_speed_travel;
network_data.direc_high_speed_travel = direc_high_speed_travel;
network_data.direc_progressions = direc_progressions;
network_data.direc_non_progressions = direc_non_progressions;


save([savedir '/' event_type '_network_data_noise_thresh_' num2str(noise_thresh) '_strength_' num2str(min_strength) '_disp_' num2str(min_disp)],'network_data')

%% Create heat maps
disp('creating heat maps...')

heat_map_with_age_sex_labels( direc_stationary_ranked, baboon_info, [-1 1], ranks_stationary )
colorbar
title([event_type ' network - stationary ' ])
print('-dpng',[fig_subdir '/' event_type '_stationary.png'])

heat_map_with_age_sex_labels( direc_low_speed_travel_ranked, baboon_info, [-1 1], ranks_low_speed_travel )
colorbar
title([event_type ' network - low-speed travel ' ])
print('-dpng',[fig_subdir '/' event_type '_low_speed_travel.png'])

heat_map_with_age_sex_labels( direc_high_speed_travel_ranked, baboon_info, [-1 1], ranks_high_speed_travel )
colorbar
title([event_type ' network - high-speed travel ' ])
print('-dpng',[fig_subdir '/' event_type '_high_speed_travel.png'])

heat_map_with_age_sex_labels( direc_progressions_ranked, baboon_info, [-1 1], ranks_progressions )
colorbar
title([event_type ' network - progressions ' ])
print('-dpng',[fig_subdir '/' event_type '_progressions.png'])

heat_map_with_age_sex_labels( direc_non_progressions_ranked, baboon_info, [-1 1], ranks_non_progressions )
colorbar
title([event_type ' network - non-progressions ' ])
print('-dpng',[fig_subdir '/' event_type '_non_progressions.png'])

heat_map_with_age_sex_labels( direc_all_ranked, baboon_info, [-1 1], ranks_all )
colorbar
title([event_type ' network - all ' ])
print('-dpng',[fig_subdir '/' event_type '_all.png'])

%% Create side-by-side comparisons

disp('creating comparison maps...')

%progressions vs non-progressions
rank_visualization( [transpose(ranks_non_progressions) transpose(ranks_progressions)], baboon_info, 'horiz' ,'non-progressions','progressions');
print('-dpng',[fig_subdir '/' event_type '_ranks_non_progressions_vs_progressions.png'])

%stationary vs high speed travel
rank_visualization( [transpose(ranks_stationary) transpose(ranks_high_speed_travel)], baboon_info, 'horiz' ,'stationary','high-speed travel');
print('-dpng',[fig_subdir '/' event_type '_ranks_stationary_vs_high_speed_travel.png'])

%all vs stationary
rank_visualization( [transpose(ranks_all) transpose(ranks_stationary)], baboon_info, 'horiz' ,'all','stationary');
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_stationary.png'])

%all vs high speed travel
rank_visualization( [transpose(ranks_all) transpose(ranks_high_speed_travel)], baboon_info, 'horiz' ,'all','high-speed travel');
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_high_speed_travel.png'])

%all vs non-progressions
rank_visualization( [transpose(ranks_all) transpose(ranks_non_progressions)], baboon_info, 'horiz' ,'all','non-progressions');
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_non_progressions.png'])

%all vs progressions
rank_visualization( [transpose(ranks_all) transpose(ranks_progressions)], baboon_info, 'horiz' ,'all','progressions');
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_progressions.png'])

%% Create regular correlation plots
%progressions vs non-progressions
plot_with_age_sex_class(transpose(ranks_non_progressions), transpose(ranks_progressions),baboon_info, 'rank in non-progressions', 'rank in progressions', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_non_progressions_vs_progressions_corr_plot.png'])

%stationary vs high speed travel
plot_with_age_sex_class(transpose(ranks_stationary), transpose(ranks_high_speed_travel),baboon_info, 'rank when stationary', 'rank during high-speed travel', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_stationary_vs_high_speed_travel_corr_plot.png'])

%stationary vs non-progressions
plot_with_age_sex_class(transpose(ranks_stationary), transpose(ranks_non_progressions),baboon_info, 'rank when stationary', 'rank in non-progressions', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_stationary_vs_non_progressions_corr_plot.png'])

%stationary vs progressions
plot_with_age_sex_class(transpose(ranks_stationary), transpose(ranks_progressions),baboon_info, 'rank when stationary', 'rank in progressions', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_stationary_vs_progressions_corr_plot.png'])

%all vs stationary
plot_with_age_sex_class(transpose(ranks_all), transpose(ranks_stationary),baboon_info, 'rank over all data', 'rank when stationary', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_stationary_corr_plot.png'])

%all vs high speed travel
plot_with_age_sex_class(transpose(ranks_all), transpose(ranks_high_speed_travel),baboon_info, 'rank over all data', 'rank during high-speed travel', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_high_speed_travel_corr_plot.png'])

%all vs non-progressions
plot_with_age_sex_class(transpose(ranks_all), transpose(ranks_non_progressions),baboon_info, 'rank over all data', 'rank during non-progressions', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_non_progressions_corr_plot.png'])

%all vs progressions
plot_with_age_sex_class(transpose(ranks_all), transpose(ranks_progressions),baboon_info, 'rank over all data', 'rank during progressions', [event_type ' network in different contexts'])
print('-dpng',[fig_subdir '/' event_type '_ranks_all_vs_progressions_corr_plot.png'])

%% Plots of pull vs anchor

%all data
plot_with_age_sex_class(transpose(pull_ranks_all), transpose(anchor_ranks_all),baboon_info, 'pull rank', 'anchor rank', 'all data')
print('-dpng',[fig_subdir '/rpull_vs_anchor_ranks_all_data_corr_plot.png'])

%stationary
plot_with_age_sex_class(transpose(pull_ranks_stationary), transpose(anchor_ranks_stationary),baboon_info, 'pull rank', 'anchor rank', 'stationary')
print('-dpng',[fig_subdir '/pull_vs_anchor_ranks_stationary_data_corr_plot.png'])

%progressions
plot_with_age_sex_class(transpose(pull_ranks_progressions), transpose(anchor_ranks_progressions),baboon_info, 'pull rank', 'anchor rank', 'progressions')
print('-dpng',[fig_subdir '/pull_vs_anchor_ranks_progressions_data_corr_plot.png'])

%non-progressions
plot_with_age_sex_class(transpose(pull_ranks_non_progressions), transpose(anchor_ranks_non_progressions),baboon_info, 'pull rank', 'anchor rank', 'non-progressions')
print('-dpng',[fig_subdir '/pull_vs_anchor_ranks_non_progressions_data_corr_plot.png'])

%high-speed travel
plot_with_age_sex_class(transpose(pull_ranks_high_speed_travel), transpose(anchor_ranks_high_speed_travel),baboon_info, 'pull rank', 'anchor rank', 'high-speed travel')
print('-dpng',[fig_subdir '/pull_vs_anchor_ranks_high_speed_travel_data_corr_plot.png'])

%low-speed travel
plot_with_age_sex_class(transpose(pull_ranks_low_speed_travel), transpose(anchor_ranks_low_speed_travel),baboon_info, 'pull rank', 'anchor rank', 'low-speed travel')
print('-dpng',[fig_subdir '/pull_vs_anchor_ranks_low_speed_travel_data_corr_plot.png'])



