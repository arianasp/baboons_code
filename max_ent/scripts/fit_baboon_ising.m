%fit Ising model to baboon data

%Inputs
state_file = '/Users/arianasp/Desktop/Baboons/output/HMM/log_sl/log_sl_HMM_data.h5';
n_days = 14; %number of days of data to use
inds_used = [1:2 4:7 9:11 15 18 19 21 22 26]; %which individuals to use in the model (the 16 who have data on all days)
weight = 0.5;
converge_thresh = 10^(-7);
max_iter = 100000;

%load data from HDF5 file
baboon_states = h5read(state_file,'/states');

%load other relevant data
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/day_start_idxs.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/baboon_info.mat')

%only include inds_used
baboon_info = baboon_info(inds_used);

fits = {};

for d = 1:n_days
    
    %pare down data to only include the individuals used, and only within the appropriate day range
    baboon_states_curr = baboon_states(inds_used,1:(day_start_idxs(d+1)-1));

    %also make states called 0 and 1 instead of 1 and 2
    baboon_states_curr = baboon_states_curr - 1;

    %get rid of any rows with NaNs
    nans = sum(isnan(baboon_states_curr),1);
    baboon_states_curr = baboon_states_curr(:,find(nans==0));
    
    d
    size(baboon_states_curr)

    [ alpha_fit, beta_fit ] = fit_ising( baboon_states, weight, converge_thresh, max_iter );
    
    fits(d).day = d;
    fits(d).alphas = alpha_fit;
    fits(d).beta = beta_fit;
    fits(d).baboon_info = baboon_info;
    fits(d).states = baboon_states_curr;

    save('/Users/arianasp/Desktop/Baboons/output/max_ent/baboon_ising_fits_by_day.mat','fits')
    disp('day completed')
end