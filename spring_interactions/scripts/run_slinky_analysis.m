%% BABOON SLINKY INTERACTIONS

%Sample workflow for getting push/pull/anchor/repel interactions between
%pairs of individuals over the course of several days, then using these
%interactions to look at patterns at the group level.

%% Required functions
%Make sure you have access to these functions (the script will not work
%otherwise)


%% Parameters (YOU CAN MODIFY THESE)

%%% ----------- Input / output parameters ----------------- %%%

%data file containing x,y data
xydatafile = '/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat';

%whether to save output data
save_data = 0;

%whether to save output figures
save_figs = 0;

%directory in which to save data (only matters if save_data == 1)
outdir = '/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data';

%directory in which to save figures (only matters if save_figs == 1)
figdir = '/Users/arianasp/Desktop/Baboons/output/push_pull/plots';

%%% ----------- Analysis-related parameters -------------- %%%

%days to include in the analysis
%typically, this should be set at as 1:14 (use first 14 days)
day_range = 1:14;

%noise threshold for picking out max / min values
%making this parameter larger will mean that larger differences in dyadic
%distance are required before they are considered events (in a sense this
%controls the spatial scale over which you want to look at interactions)
noise_thresh = 10;

%which type of interaction to focus on (push, pull, anchor, repel)
event_type = 'pull';

%which type of threshold to use (multiplicative or additive) 'mult' or
%'add'
thresh_type = 'mult';

%minimum strength of interaction (multiplicative strength) needed to
%count it as an event
strength_thresh = 0.2;

%minimum disparity of interaction (multiplicative disparity) needed to
%count it as an event
disp_thresh = 0.2;

%% Load required data

disp('loading data...')

load(xydatafile)


%% Get interactions

disp('getting interactions...')

N = size(xs,1); %number of individuals

interactions_all = cell(N,N); %initialize cell array to hold data

%get interactions between every pair of individuals and store them
for a = 1:N
    for b = (a+1):N
    	[ interactions ] = dyadic_interactions( xs, ys, day_start_idxs, day_range, a, b, noise_thresh );
        interactions_all{a,b} = interactions;
	end
end

%% Construct interaction matrices
%event_mat(i,j) is the number of interactions in which j lead i.

%direc_mat is the overall directionality of the relationship between 
%i and j.  This is a number that ranges from -1 (i leads more often) 
%through 0 (neither leads more often) to 1 (j leads more often). It is
%defined as: direc_mat(i,j) = [event_mat(i,j) - 
%               event_mat(j,i)]/[event_mat(i,j)+ event_mat(j,i)]

disp('constructing event matrix...')

[ event_mat, direc_mat ] = event_count_matrix( interactions_all, event_type, thresh_type, strength_thresh, disp_thresh );

%% Rank directionality matrix by mean directionality for each individual
%Rank the directionality matrix (direc_mat) from high to low. "Leaders"
%will be displayed at the bottom row


disp('ranking event matrix...')
[ new_mat, ranks ] = rank_mat( direc_mat, 'mean', 'descend', 'col');
for i = 1:N
    new_mat(i,i) = 0;
end


%% Create figure
disp('creating figure...')

heat_map_with_age_sex_labels( new_mat, baboon_info, [-1 1], ranks )




