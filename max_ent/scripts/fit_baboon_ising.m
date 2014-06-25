%fit Ising model to baboon data

%Inputs
state_file = '/Users/arianasp/Desktop/Baboons/output/HMM/log_sl/log_sl_HMM_data.h5';
n_days = 14; %number of days of data to use
inds_used = [1:7 9:11 15 18 19 21 22 26]; %which individuals to use in the model (the 16 who have data on all days)
weight = 0.5;
converge_thresh = 10^(-7);
max_iter = 100000;

%load data from HDF5 file
baboon_states = h5read(state_file,'/states');

%load other relevant data
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/day_start_idxs.mat')
load('/Users/arianasp/Desktop/Baboons/data/matlab_raw/baboon_info.mat')

%pare down data to only include the individuals used, and only within the appropriate day range
baboon_states = baboon_states(inds_used,1:(day_start_idxs(n_days+1)-1));
baboon_info = baboon_info(inds_used);

%also make states called 0 and 1 instead of 1 and 2
baboon_states = baboon_states - 1;

%get rid of any rows with NaNs
nans = sum(isnan(baboon_states),1);
baboon_states = baboon_states(:,find(nans==0));

size(baboon_states)

[ alpha_fit, beta_fit ] = fit_ising( data, weight, converge_thresh, max_iter );