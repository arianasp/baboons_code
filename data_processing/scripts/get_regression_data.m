%create a dataframe containing information for a regression, which will be
%entered into R
%columns: individual ID, sex (0 = female, 1 = male), age (0 = adult, = 1
%subadult, 2 = juvenile), weight, day, timestep, dist from centroid, dist
%in front of centroid, normalized dist from centroid, normalized dist in
%front of centroid

clear all;

downsample_rate = 60;

max_day = 14;

disp('loading data')


load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/terrain_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/path_direcs_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/pos_rel_to_centroid_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/sleep_sites_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/pos_rel_to_centroid_ranked_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/troop_level_metrics/troop_pers_vel_step10_level1.mat','troop_pers_vels')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/speeds_step_10_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/mahal_dists_level1.mat','mahal_dists')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics/mean_dist_to_others_level1.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_processed/troop_level_metrics/order_params_level1.mat')


disp('done loading data:')
whos

lateral_dists = abs(xs_rel);

%fix weight of individual 15 to 8 (NOTE: IF USING THIS FOR DIFFERENT DATA SET,
%NEED TO TAKE THIS HARD-CODING OUT!)
baboon_info(15).weight = {'8'};

%filter data
%put nans in where there are fewer than 10 tracked individuals and less
%than 100 m from sleep site 
n_inds_tracked = sum(not(isnan(xs)),1);
times_to_include = find(n_inds_tracked >= 10);
times_to_include2 = find(troop_dist_from_sleep_site >= 100);
times_to_include = intersect(times_to_include,times_to_include2);
times_to_include = intersect(times_to_include,1:(day_start_idxs(max_day+1)-1));

%get basic information about how big matrix needs to be
n_days = length(day_start_idxs);
n_inds = length(baboon_info);
n_times = size(xs,2);
day_start_idxs = [day_start_idxs n_times+1];

%create matrix to hold all data
n_rows = length(times_to_include)*n_inds;
n_cols = 20;
data = nan(n_rows,n_cols);

%fill matrix with data
idx = 1;
t_prev = 1;
times_to_include(end)
for time_idx = 1:length(times_to_include)
    t = times_to_include(time_idx);
    d = find(day_start_idxs <= t,1,'last');
    if (t - t_prev) >= downsample_rate
        t_prev = t;
        t
        for i = 1:n_inds
            if not(isnan(dist_from_centroid(i,t)))
                data(idx,1) = d; %1. day
                data(idx,2) = t; %2. time step idx
                data(idx,3) = baboon_info(i).animal_id; %3. animal id
                data(idx,4) = strcmp(baboon_info(i).sex,'M'); %4. sex - 0 for female, 1 for male
                if strcmp(baboon_info(i).age,'A') %5. age - 0 for adult, 1 for subadult, 2 for juv
                    data(idx,5) = 0;
                elseif strcmp(baboon_info(i).age,'SA')
                    data(idx,5) = 1;
                else
                    data(idx,5) = 2;
                end
                data(idx,6) = str2num(baboon_info(i).weight{1}); %6. weight
                data(idx,7) = troop_path_direcs301(t); %7. troop path direc
                data(idx,8) = troop_pers_vels(t); %8. troop pers vel
                data(idx,9) = troop_speeds(t); %9. troop speed
                data(idx,10) = dist_from_centroid(i,t); %10. ind dist from centroid
                data(idx,11) = ys_rel(i,t); %11. distance in front of centroid
                data(idx,12) = mahal_dists(i,t); %14. mahalanobis distance from centroid
                data(idx,13) = terrain(i,t); %individual terrain
                data(idx,14) = troop_terrain(t); %troop most common terrain
                data(idx,15) = speeds(i,t); %individual speed
                data(idx,16) = mean_dist_to_others(i,t); %mean distance to other individuals
                data(idx,17) = lateral_dists(i,t); %lateral distance from centroid i.e. abs(rel_xs) 
                data(idx,18) = order_params.eccentricity(t); %eccentricity
                data(idx,19) = order_params.tilt(t); %tilt
                data(idx,20) = order_params.polarization(t); %polarization
                idx = idx + 1;
            end
        end
        t
    end
end

disp('done gathering data...')

%remove excess NaNs
data = data(1:idx,:);

disp('saving data...')

%save to HDF5
h5create('/Users/arianasp/Desktop/Dropbox/baboons_shared/positioning_analysis/data/regression_data_days_1-14_min_inds_10_sleep_dist_100_downsample_60.h5','/data',size(data))
h5write('/Users/arianasp/Desktop/Dropbox/baboons_shared/positioning_analysis/data/regression_data_days_1-14_min_inds_10_sleep_dist_100_downsample_60.h5','/data',data)

    
    
    
    
    
    
    
    
    
    
    
    
    
    