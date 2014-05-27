function [ times, strengths, disparities] = interaction_stats( interactions_all, event_type, make_plots, a, b )
%Gets lengths (in time), strengths, and disparities for all interaction
%events of a certain type (event_type). If a and b are specified, it will
%get interactions only between individuals a and b. If make_plots == 1, it
%will output a plot giving the distributions of these three measures.
%INPUTS:
%   interactions_all: [N x N cell array] containing structs with
%       information about interaction events (see run_slinky_analysis for
%       how to generate this cell array)
%   event_type: [string] telling which events to pull out ('pull,'anchor',
%       'repel','push')
%   make_plots: [bool] whether to produce plots or not (defaults to 0)
%   a,b: individuals to get data about - a is the leader, b the follower
%       (if this is unspecified, interactions will be aggregated over all 
%       pairs of individuals
%OUTPUTS:
%   times: [n_ints x length(event_types) matrix] how long each interaction took
%   strengths: [n_ints x length(event_types) matrix] strength of each
%       interaction
%   disparities: [n_ints x length(event_types) matrix] disparity of each
%       interaction

if ~exist('make_plots')
    make_plots = 0;
end

if exist('a') && exist('b')
    all_inds = 0;
else
    all_inds = 1;
end


times = []; strengths = []; disparities = [];
if all_inds
    for i = 1:size(interactions_all,1)
        for j = (i+1):size(interactions_all,2)
            interactions = interactions_all{i,j};
            for k = 1:length(interactions)
                if strcmp(interactions(k).type,event_type) 
                    times = [times interactions(k).time_idxs(3) - interactions(k).time_idxs(1)];
                    strengths = [strengths interactions(k).strength_mult];
                    disparities = [disparities interactions(k).disp_mult];
                end
            end
        end
    end
else
    interactions = interactions_all{min(a,b),max(a,b)};
    for k = 1:length(interactions)
        if interactions(k).leader == a
            if strcmp(interactions(k).type,event_type)
                times = [times interactions(k).time_idxs(3) - interactions(k).time_idxs(1)];
                strengths = [strengths interactions(k).strength_mult];
                disparities = [disparities interactions(k).disp_mult];
            end
        end
    end
end

if make_plots
    figure
    stren_hist = histc(strengths,0:.05:1);
    disp_hist = histc(disparities,0:.05:1);
    time_hist = histc(log10(times),0:.05:4);
    stren_hist = stren_hist / sum(stren_hist);
    disp_hist = disp_hist / sum(disp_hist);
    time_hist = time_hist / sum(time_hist);
    subplot(1,3,1)
    bar(0:.05:4,time_hist,'histc')
    xlabel('log10(time) (sec)')
    ylabel('probability')
    title('TIME')
    xlim([0 4])
    subplot(1,3,2)
    bar(0:.05:1,stren_hist,'histc')
    xlabel('strength')
    title('STRENGTH')
    xlim([0 1])
    subplot(1,3,3)
    bar(0:.05:1,disp_hist,'histc')
    xlabel('disparity')
    title('DISPARITY')
    xlim([0 1])
end
    
    

    



end

