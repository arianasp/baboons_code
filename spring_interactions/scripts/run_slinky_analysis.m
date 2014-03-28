%% BABOON SLINKY INTERACTIONS

%Sample workflow for getting push/pull/anchor/repel interactions between
%pairs of individuals over the course of several days, then using these
%interactions to look at patterns at the group level.

%% Required functions
%Make sure you have access to these functions (the script will not work
%otherwise)


%% Parameters (YOU CAN MODIFY THESE)

%%% ----------- Input / output parameters ----------------- %%%

%whether to load data from file or generate new data
load_data = 1;

%whether to save output data
save_data = 0;

%whether to save output figures
save_figs = 0;

%data file containing x,y data
xydatafile = '/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat';

%directory in which to save data (only matters if save_data == 1)
outdir = '/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data';

%directory in which to save figures (only matters if save_figs == 1)
figdir = '/Users/arianasp/Desktop/Baboons/output/push_pull/plots';

datadir = '/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data';

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
strength_thresh = 0.1;

%minimum disparity of interaction (multiplicative disparity) needed to
%count it as an event
disp_thresh = 0.1;

%% Load required data

disp('loading data...')

if load_data
    %load pre-computed interactions data
    load([datadir '/dyad_interactions_noise_thresh_' num2str(noise_thresh) '.mat'])
    load(xydatafile)
    N = size(interactions_all,1); %number of individuals
else
    %load raw x-y coordinate data to compute interactions now
    load(xydatafile)
    N = size(xs,1);
end


%% Get interactions if needed

if not(load_data)
    disp('getting interactions...')

    interactions_all = cell(N,N); %initialize cell array to hold data

    %get interactions between every pair of individuals and store them
    for a = 1:N
        for b = (a+1):N
            [ interactions ] = dyadic_interactions( xs, ys, day_start_idxs, day_range, a, b, noise_thresh );
            interactions_all{a,b} = interactions;
        end
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

%create a heat map
heat_map_with_age_sex_labels( new_mat, baboon_info, [-1 1], ranks )
colorbar
title([event_type ' directionality, noise thresh = ', num2str(noise_thresh), ', min mult. strength = ' num2str(strength_thresh) ', min disparity = ' num2str(disp_thresh)])

if save_figs
    print('-dpng',[figdir '/' event_type '_direc_mat_noise_thresh_' num2str(noise_thresh) '_min_strength_' num2str(strength_thresh) '_min_disp_' num2str(disp_thresh) '.png'])
end

if save_data
    save([outdir '/dyad_interactions_noise_thresh_' num2str(noise_thresh) '.mat'])
end


