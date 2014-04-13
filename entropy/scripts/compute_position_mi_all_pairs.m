%get spatial mutual information for all pairs of individuals

%% Parameters
nbins = [20 20];
mins = [-100 -100];
maxes = [100 100];
dir = '/Users/samiam/Desktop/Baboons';
ndays = 14;

paras.nbins = nbins;
paras.mins = mins;
paras.maxes = maxes;
paras.ndays = ndays;

%% Load data
load([dir '/data/matlab_processed/individual_level_metrics/pos_rel_to_centroid_level1.mat'])
load([dir '/data/matlab_raw/day_start_idxs.mat'])

xs = xs_rel(:,1:(day_start_idxs(ndays+1)-1));
ys = ys_rel(:,1:(day_start_idxs(ndays+1)-1));

N = size(xs_rel,1);

%% Compute MI for all pairs
MI_data = cell(N,N);
MI_mat = nan(N,N);
for i = 1:N
    i
    for j = i:N
        X = [xs(i,:)' ys(i,:)'];
        Y = [xs(j,:)' ys(j,:)'];
        dataij = mutual_information(X,Y,mins,maxes,nbins);
        MI_data{i,j} = dataij;
        MI_data{j,i} = dataij;
        MI_mat(i,j) = dataij.MI;
        MI_mat(j,i) = dataij.MI;
    end
end

save([dir '/position_mi/position_mi_days1-' num2str(ndays) '.mat'],'MI_data','MI_mat','paras')
        