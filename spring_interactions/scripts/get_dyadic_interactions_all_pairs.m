%get dyadic interactions for all pairs and construct cell arrays of structs

noise_thresh_vals = [3 5 8 12 15 20 25 30 35 40 50];
day_range = 1:14;
outdir = '/Users/arianasp/Desktop/Baboons/push_pull/interactions_data';
datafile = '/Users/arianasp/Desktop/Baboons/data/matlab_raw/xy_level1.mat';

load(datafile)
N = size(xs,1);

for noise_thresh = noise_thresh_vals
    interactions_all = cell(N,N);
    for a = 1:N
        a
        for b = (a+1):N
            b
            [ interactions ] = dyadic_interactions( xs, ys, day_start_idxs, day_range, a, b, noise_thresh );
            interactions_all{a,b} = interactions;
        end
    end
    save([outdir '/dyad_interactions_noise_thresh_' num2str(noise_thresh) '.mat'],'interactions_all')
end