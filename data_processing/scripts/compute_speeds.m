%Compute and save speeds from positions data.

%% Parameters

%time step to use when computing the speeds
%if this is equal to 1, the speed is found using the difference between
%positions in consecutive time steps. If this is greater than 1, the speed
%is found using the difference between positions that are t_step seconds
%apart.
t_step = 3;

%filename for the file with xy data
xy_file = ['/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat'];

%filename for the file with centroid over time data
centr_file = ['/Users/arianasp/Desktop/Baboons/data/matlab_processed/troop_level_metrics/centroid_min_inds_10_level1.mat'];

%directory in which to save the output file
savedir = ['/Users/arianasp/Desktop/Baboons/data/matlab_processed/individual_level_metrics'];

%% Load files

disp('loading files...')

load(xy_file)
load(centr_file)

%% Extract appropriate data

disp('extracting data...')

%number of days
n_days = length(day_start_idxs);

%add an ending time to day_start_idxs
day_start_idxs = [day_start_idxs size(xs,2)+1];

%position data (regular and shifted up by t_step)
%note the precision change to "double" - this is because we are subtracting
%small numbers to get speeds
xs = double(xs);
xs_shift = xs(:,(1+t_step):end);
xs_curr = xs(:,1:(end-t_step));
ys = double(ys);
ys_shift = ys(:,(1+t_step):end);
ys_curr = ys(:,1:(end-t_step));

troopxs = double(centr_over_time(:,1));
troopxs_shift = troopxs((1+t_step):end);
troopxs_curr = troopxs(1:(end-t_step));
troopys = double(centr_over_time(:,2));
troopys_shift = troopys((1+t_step):end);
troopys_curr = troopys(1:(end-t_step));

%% Compute speeds

disp('computing speeds...')

%x and y velocities for individuals
dxs = (xs_shift - xs_curr)/t_step; 
dys = (ys_shift - ys_curr)/t_step;

%x and y velocities for troop centroid
troopdxs = (troopxs_shift - troopxs_curr)/t_step; 
troopdys = (troopys_shift - troopys_curr)/t_step;

%number of individuals tracked over time
n_inds = sum(not(isnan(xs)),1); 
n_inds_shift = n_inds((1+t_step):end);
n_inds_curr = n_inds(1:(end-t_step));

%change in number of tracked individuals (need to make troop speeds 
%adjacent to these NaN)
d_n_inds = n_inds_shift - n_inds_curr; 

%get rid of speeds that cross days (replace with NaNs)
for d = 1:n_days
    last = day_start_idxs(d+1) - 1;
    dxs(:,(last-t_step+1):last) = NaN;
    dys(:,(last-t_step+1):last) = NaN;
    troopdxs((last-t_step+1):last) = NaN;
    troopdys((last-t_step+1):last) = NaN;
end

%get rid of troop-level speeds where the number of tracked individuals
%changed (replace with NaNs)
changing_inds = find(abs(d_n_inds)>0);
for i = 1:length(changing_inds)
    idx = changing_inds(i);
    troopdxs(idx) = NaN;
    troopdys(idx) = NaN;
end

%compute speeds for individual and troop level
troop_speeds = sqrt(troopdxs.^2 + troopdys.^2);
speeds = sqrt(dxs.^2 + dys.^2);

%% Save data

disp('saving data...')

savename = [savedir '/speeds_step_' num2str(t_step) '_level1.mat'];
save(savename,'speeds','troop_speeds');

disp('done')






