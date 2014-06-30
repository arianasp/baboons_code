%get a more reliable estimate of the entropy by extrapolating from
%estimates based on subsamples

%Inputs
state_file = '/home/arianasp/baboons/output/HMM/log_sl_HMM_data.h5';
datadir = '/home/arianasp/baboons/data/matlab_raw';
outdir = '/home/arianasp/baboons/output/max_ent';
n_days = 14; %number of days of data to use
inds_used = [1:2 4:7 9:11 15 18 19 21 22 26]; %which individuals to use in the model (the 16 who have data on all days)
data_fracs = linspace(1,10,19);
reps = 100;


%READ IN AND PRE-PROCESS DATA

%load data from HDF5 file
baboon_states = h5read(state_file,'/states');

%load other relevant data
load([datadir '/day_start_idxs.mat'])
load([datadir '/baboon_info.mat'])

%only include inds_used
baboon_info = baboon_info(inds_used);

%pare down data to only include the individuals used, and only within the appropriate day range
baboon_states_curr = baboon_states(inds_used,1:(day_start_idxs(n_days+1)-1));

%also make states called 0 and 1 instead of 1 and 2
baboon_states_curr = baboon_states_curr - 1;

%get rid of any rows with NaNs
nans = sum(isnan(baboon_states_curr),1);
baboon_states_curr = baboon_states_curr(:,find(nans==0));

%SUBSAMPLE AND ESTIMATE ENTROPY
[ H_ests ] = entropy_vs_sample_size( baboon_states_curr, data_fracs, reps )

save([outdir '/entropy_vs_sample_size_100reps.mat'])