function [ event_mat, direc_mat ] = event_count_matrix( interactions_all, event_type, measure_type, strength_thresh, disp_thresh, times_to_include )
%Constructs a matrix mat, where mat(i,j) = # of events of type "event_type"
%where j lead i (i.e. either j pulled i, j pushed i, j anchored i, or j
%repelled i).
%INPUTS:
%   interactions_all: [NxN cell array] containing structs with all 
%       interactions between each pair of individuals 
%       (output of get_dyadic_interactions_all_pairs)
%   event_type: [string] specifying which type of event (can be
%       "push","pull","anchor" or "repel"
%   measure_type: [string] specifying whether to use the multiplicative 
%       ("mult")or additive ("add") version of the disparity and strength 
%       measures
%   strength_thresh: [number] minimum strength needed to include an event 
%       in the count
%   disp_thresh: [number] minimum disparity needed to include an event 
%       in the count
%   times_to_include: [vector] of time indexes that should be included when
%       counting events. This is so that you can create a matrix that
%       represents the number of interactions across a subset of the data.
%       By default, times_to_include includes all times.
%OUTPUTS:
%   event_mat: [NxN matrix] where mat(i,j) = number of times that j lead i
%   direc_mat: [NxN matrix] where mat(i,j) = directionality of the
%       relationship between i and j. This is a number that ranges from -1
%       (i leads more often) through 0 (neither leads more often) to 1 (j
%       leads more often). It is defined as:
%       direc_mat(i,j) = [event_mat(i,j) - event_mat(j,i)]/[event_mat(i,j)
%       + event_mat(j,i)]

%get number of individuals
N = size(interactions_all,1);

%initialize matrix to hold results
mat = zeros(N,N);

if exist('times_to_include')
    subset_time = 1;
else
    subset_time = 0;
end

%for each pair of individuals
for i = 1:N
    i
    for j = (i+1):N
        %get interactions between those individuals
        interactions = interactions_all{i,j};
        for k = 1:length(interactions)
            
            if strcmp(interactions(k).type,event_type)
                %get leader and follower
                lead = interactions(k).leader;
                foll = interactions(k).follower;
            
                %get disparity and strength
                if strcmp(measure_type,'mult')
                    disp = interactions(k).disp_mult;
                    strength = interactions(k).strength_mult;
                elseif strcmp(measure_type,'add')
                    disp = interactions(k).disp_add;
                    strength = interactions(k).strength_add;
                else
                    error('must specify measure_type as "mult" or "add"')
                end
                
                %if disparity and strength are high enough, include the
                %event in the count (add 1 to the corresponding matrix
                %cell)
                if disp >= disp_thresh && strength >= strength_thresh
                    if not(subset_time)
                        mat(foll,lead) = mat(foll,lead) + 1;
                    else
                       if any(times_to_include == interactions(k).time_idxs(2))
                           mat(foll,lead) = mat(foll,lead) + 1;
                       end
                    end
                end
            end
        end
    end
end

%assign the matrix of event counts to event_mat
event_mat = mat;

%initialize a matrix to hold the directionality information
direc_mat = zeros(N,N);

%create the directionality matrix
for i = 1:N
    for j = 1:N
        if i~=j
            direc_mat(i,j) = (event_mat(i,j) - event_mat(j,i))/(event_mat(i,j)+event_mat(j,i));
        else
            direc_mat(i,j) = NaN;
        end
    end
end
                
            
            


end

