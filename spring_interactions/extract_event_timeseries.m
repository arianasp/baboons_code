function [ interactions_over_time ] = extract_event_timeseries( interactions_all, event_type, t_start, t_end, min_strength, min_disp )
%Get a timeseries of events that occur between pairs of individuals, and
%store everything in an NxNxT array.
%INPUTS:
%   interactions_all: [cell array] containing interaction information 
%       (contained in structs) for each pair of individuals 
%       (output of get_dyadic_interactions_all_pairs script)
%   event_type: [string] representing what type of events to pull out - can
%       be 'push','pull','anchor', or 'repel'
%   t_start: [number] the first time index to include
%   t_end: [number] the last time index to include
%   min_strength: [number] minimum strength required to count an event
%       (defaults to 0)
%   min_disp: [number] minimum disaprity required to count an event
%       (defaults to 0)
%OUTPUTS:
%   interactions_over_time: [NxNx(t_end - t_start + 1) array] where
%       interactions_over_time(i,j,t) = 1 if j is leading i at time t and 0
%       otherwise

%get number of individuals
N = size(interactions_all,1);

%check that the data is in the correct form (NxN cell array)
if size(interactions_all,2) ~= N
    error('interactions_all must be a square cell array')
end

%check that the event type is one of the supported ones
if not(strcmp(event_type,'push') || strcmp(event_type,'pull') || ...
        strcmp(event_type,'anchor') || strcmp(event_type,'repel'))
    error('event_type must be "push","pull","anchor", or "repel"')
end

%default min_strength = 0
if ~exist('min_strength')
    min_strength = 0;
end

%default min_disp = 0
if ~exist('min_disp')
    min_disp = 0;
end

%create an array to hold timeseries data
interactions_over_time = zeros(N,N,(t_end - t_start + 1));

%get interactions and put them in the matrix
for i = 1:N
    for j = (i+1):N
        
        %get interactions for each pair
        interactions = interactions_all{i,j};
        
        %for each interaction
        for k = 1:length(interactions)
            
            t0 = interactions(k).time_idxs(1);
            tf = interactions(k).time_idxs(3);
            
            %if within the time range
            if t0 >= t_start || tf <= t_end
                
            
                %get its type
                curr_type = interactions(k).type;
            
                %count it if it's the correct type and passes the thresholds
                if strcmp(curr_type,event_type)
                    strength = interactions(k).strength_mult;
                    disp = interactions(k).disp_mult;
                    if strength >= min_strength && disp >= min_disp
                        t_event_begin = max(t0,t_start);
                        t_event_end = min(tf,t_end);
                        start_idx = t_event_begin - t_start + 1;
                        end_idx = t_event_end - t_start + 1;
                        interactions_over_time(i,j,start_idx:end_idx) = 1;
                    end
                end
            end
        end
    end
end

        


end

