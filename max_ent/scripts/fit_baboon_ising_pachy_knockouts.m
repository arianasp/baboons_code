%fit Ising model to baboon data

%Inputs
state_file = '/home/arianasp/baboons/output/HMM/log_sl_HMM_data.h5';
datadir = '/home/arianasp/baboons/data/matlab_raw';
outdir = '/home/arianasp/baboons/output/max_ent';
n_days = 14; %number of days of data to use
inds_used = [1:2 4:7 9:11 15 18 19 21 22 26]; %which individuals to use in the model (the 16 who have data on all days)
weight = 0.5;
converge_thresh = 10^(-7);
max_iter = 100000;

%load data from HDF5 file
baboon_states = h5read(state_file,'/states');

%load other relevant data
load([datadir '/day_start_idxs.mat'])
load([datadir '/baboon_info.mat'])

%only include inds_used
baboon_info = baboon_info(inds_used);

fits = {};

for i = 1:length(inds_used)
    
    inds_used_curr = inds_used;
    inds_used_curr(i) = []; %knock out individual inds_used(i)
    
    baboon_info_used = baboon_info;
    baboon_info_used(i) = [];
    baboon_info_knockout = baboon_info_used(i);
    
    %pare down data to only include the individuals used, and only within the appropriate day range
    baboon_states_curr = baboon_states(inds_used_curr,1:(day_start_idxs(n_days+1)-1));

    %also make states called 0 and 1 instead of 1 and 2
    baboon_states_curr = baboon_states_curr - 1;

    %get rid of any rows with NaNs
    nans = sum(isnan(baboon_states_curr),1);
    baboon_states_curr = baboon_states_curr(:,find(nans==0));
    
    inds_used(i)
    inds_used_curr

    [ alpha_fit, beta_fit ] = fit_ising( baboon_states_curr, weight, converge_thresh, max_iter );
    
    fits(i).inds_used = inds_used_curr;
    fits(i).knocked_out = inds_used(i);
    fits(i).alphas = alpha_fit;
    fits(i).beta = beta_fit;
    fits(i).baboon_info_used = baboon_info_used;
    fits(i).baboon_info_knockout = baboon_info_knockout;
    fits(i).states = baboon_states_curr;

    save([outdir '/baboon_ising_fits_single_knockouts.mat','fits'])
    disp('ind completed')
end