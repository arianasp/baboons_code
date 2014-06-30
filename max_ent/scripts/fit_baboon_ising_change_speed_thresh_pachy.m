%fit Ising model to baboon data

%Inputs
n_days = 14; %number of days of data to use
inds_used = [1:2 4:7 9:11 15 18 19 21 22 26]; %which individuals to use in the model (the 16 who have data on all days)
weight = 0.5;
converge_thresh = 10^(-7);
max_iter = 100000;
threshes = [-4 -3.5 -3 -2.5 -2:.2:0 0.5];

%load data from speeds file
load('/home/arianasp/baboons/data/matlab_processed/individual_level_metrics/speeds_step_10_level1.mat')
logspeeds = log(speeds);

%load other relevant data
load('/home/arianasp/baboons/data/matlab_raw/day_start_idxs.mat')
load('/home/arianasp/baboons/data/matlab_raw/baboon_info.mat')

%only include inds_used
baboon_info = baboon_info(inds_used);

fits = {};

for i = 1:length(threshes)
    thresh = threshes(i);
    
    baboon_states = logspeeds > thresh;
    
    %pare down data to only include the individuals used, and only within the appropriate day range
    baboon_states_curr = baboon_states(inds_used,1:(day_start_idxs(n_days+1)-1));

    %get rid of any rows with NaNs
    nans = sum(isnan(speeds(inds_used,:)),1);
    baboon_states_curr = baboon_states_curr(:,nans==0);
    
    size(baboon_states_curr)

    [ alpha_fit, beta_fit ] = fit_ising( baboon_states_curr, weight, converge_thresh, max_iter );
    
    fits(i).log_speed_thresh = thresh;
    fits(i).alphas = alpha_fit;
    fits(i).beta = beta_fit;
    fits(i).inds_used = inds_used;
    fits(i).n_days = n_days;
    fits(i).baboon_info = baboon_info;
    fits(i).states = baboon_states_curr;

    save('/home/arianasp/baboons/output/max_ent/baboon_ising_fits_by_speed_thresh3.mat','fits')
    disp('thresh completed')
end