function [ ] = interactions_matlab_to_csv( dir, filename  )
%Takes in the filename of an interactions cell array (in Matlab form) 
%and converts it into a csv for import into R. Each row of the csv file 
%contains one interaction, and the columns are as follows:
%   day: day that the interaction occurred
%   noise_thresh: noise threshold used to pull out the interactions
%   t1: start time of the interaction
%   t2: middle time of the interaction
%   t3: end time of the interaction
%   strength_add: additive strength of the interaction
%   strength_mult: multiplicative strength of the interaction
%   disp_add: additive disparity of the interaction
%   disp_mult: multiplicative disparity of the interaction
%   type: type of interaction
%   leader: index of the leading individual
%   follower: index of the following individual
%INPUTS:
%   dir: directory where the input file is
%       e.g.'/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data'
%   filename: name of the file (without extension)
%       e.g. 'dyad_interactions_noise_thresh_10'
%OUTPUTS:
%   none, but saves a file in the same directory (dir), with the same
%       filename but the extension csv, in the format specified above

%load data
disp('loading data...')
load([dir '/' filename '.mat'])

if ~exist('interactions_all')
    error('no variable called interactions_all exists - check the input data file')
end

%number of individuals
N = size(interactions_all,1);

%name of the output file
outfile = [dir '/' filename '.csv'];

%check if the file already exists
disp('checking if output file exists already...')
if exist(outfile,'file')==2
    error('output file already exists - delete the file and try again')
end

%create and open the output file
disp('creating output file...')
fid = fopen(outfile,'wt');

%create a header
fprintf(fid,'%c','day,noise.thresh,t1,t2,t3,strength.add,strength.mult,disp.add,disp.mult,type,leader,follower')
fprintf(fid,'\n')

disp('writing data to file...')
for i = 1:N
    i
    for j = (i+1):N
        j
        interactions = interactions_all{i,j};
        for k = 1:length(interactions)
            fprintf(fid,'%c',num2str(interactions(k).day));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).noise_thresh));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).time_idxs(1)));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).time_idxs(2)));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).time_idxs(3)));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).strength_add));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).strength_mult));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).disp_add));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).disp_mult));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',interactions(k).type);
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).leader));
            fprintf(fid,'%c',',');
            fprintf(fid,'%c',num2str(interactions(k).follower));
            fprintf(fid,'\n');
        end
    end
end

disp('closing file...')
fclose(fid);
disp('done')
        



end

