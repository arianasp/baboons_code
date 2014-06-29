%fit Ising model to baboon data

%Inputs
n_days = 14; %number of days of data to use
inds_used = [1:2 4:7 9:11 15 18 19 21 22 26]; %which individuals to use in the model (the 16 who have data on all days)
weight = 0.5;
converge_thresh = 10^(-7);
max_iter = 100000;
threshes = [-3 -3 -1.9 -1.8 -1.7 -1.6 -1.5 -1.4 -1.3 -1.2 -1.1 -1 -.9 -.8 -.7 -.6 -.5 -.4 -.3 -.2 -.1 0 0.1 0.2 0.3 0.4 0.5];

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
    
    baboon_states = speeds > thresh;
    
    %pare down data to only include the individuals used, and only within the appropriate day range
    baboon_states_curr = baboon_states(inds_used,1:(day_start_idxs(n_days+1)-1));

    %get rid of any rows with NaNs
    nans = sum(isnan(baboon_states_curr),1);
    baboon_states_curr = baboon_states_curr(:,find(nans==0));
    
    size(baboon_states_curr)

    [ alpha_fit, beta_fit ] = fit_ising( baboon_states_curr, weight, converge_thresh, max_iter );
    
    fits(i).log_speed_thresh = thresh;
    fits(i).alphas = alpha_fit;
    fits(i).beta = beta_fit;
    fits(i).inds_used = inds_used;
    fits(i).n_days = n_days;
    fits(i).baboon_info = baboon_info;
    fits(i).states = baboon_states_curr;

    save('/home/arianasp/baboons/output/max_ent/baboon_ising_fits_by_speed_thresh.mat','fits')
    disp('thresh completed')
end