function [  ] = interactions_over_time_animation( interactions_all, baboon_info, event_type, t_start, t_end, t_step, min_strength, min_disp, imgdir )
%Outputs images for the construction of an animation of interactions over
%time. 
%INPUTS:
%   interactions_all: [cell array] containing interaction information 
%       (contained in structs) for each pair of individuals 
%       (output of get_dyadic_interactions_all_pairs script)
%   event_type: [string] representing what type of events to pull out - can
%       be 'push','pull','anchor', or 'repel'
%   t_start: [number] the first time index to include
%   t_end: [number] the last time index to include
%   t_step: [number] time step (e.g. if 5, show an image for every 5 sec)
%   min_strength: [number] minimum strength required to count an event
%       (defaults to 0)
%   min_disp: [number] minimum disaprity required to count an event
%       (defaults to 0)
%   imgdir: directory in which to save the images
%OUTPUTS:
%   png images of an NxN matrix for every time step, where black squares
%   indicate interactions

[ interactions_over_time ] = extract_event_timeseries( interactions_all, event_type, t_start, t_end, min_strength, min_disp );

T = size(interactions_over_time,3);
N = size(interactions_over_time,1);

cd(imgdir)
dir = [imgdir '/' event_type 't_' num2str(t_start) '-' num2str(t_end) 'step_' num2str(t_step) '_min_stren_' num2str(min_strength) '_min_disp_' num2str(min_disp)];
mkdir(dir);

for t = 1:t_step:T
    heat_map_with_age_sex_labels( interactions_over_time(:,:,t), baboon_info, [-1 1], 1:N );
    xlabel(['t = ' num2str(t_start + t - 1)])
    title([event_type ' events'])
    print('-dpng',[dir '/' num2str(t) '.png'])
    close(gcf)
end
    


end

