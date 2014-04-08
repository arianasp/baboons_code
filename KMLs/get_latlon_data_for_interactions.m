%generate csv files with lat/lon data for each interaction

noise_thresh = 10;
interaction_type = 'pull';
max_ints_per_pair = 10;

gps_dir = '/Users/arianasp/Desktop/Baboons/data/matlab_raw';
interactions_dir = '/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data';
out_dir = ['/Users/arianasp/Desktop/Baboons/google_earth/slinky_interactions/noise_thresh_' num2str(noise_thresh) '/latlons'];

mkdir(out_dir)

disp('loading data...')
%load gps data
load([gps_dir '/gps_processed_level1.mat'])

%load interactions data
load([interactions_dir '/' 'dyad_interactions_noise_thresh_' num2str(noise_thresh) '.mat'])

%number of individuals
N = size(interactions_all,1);

N=2;

%file to hold data about each interaction
listfid = fopen([out_dir '/interactions_list.csv'],'wt')
fprintf(listfid,'FILE,T1,T2,LEADER,FOLLOWER\n')

disp('generating csv files...')
format long
for i = 1:N
	i
	for j = (i+1):N
		j
		interactions = interactions_all{i,j};
		for k = 1:min(max_ints_per_pair, length(interactions))
			k
			t1 = interactions(k).time_idxs(1);
			t2 = interactions(k).time_idxs(3);
			ts = t1:t2;
			lead = interactions(k).leader;
			foll = interactions(k).follower;
			strength = interactions(k).strength_mult;
			disp = interactions(k).disp_mult;
			filename = [interaction_type '_lead_' num2str(lead) '_foll_' num2str(foll) '_noise_' num2str(noise_thresh) '_strength_' num2str(strength) '_disp_' num2str(disp) '_t_' num2str(t1) '-' num2str(t2)];
			fprintf(listfid,[filename ',' num2str(t1) ',' num2str(t2) ',' num2str(lead) ',' num2str(foll) '\n'])
			fid = fopen([out_dir '/' filename '.csv'],'wt');
			fprintf(fid,'ID,TIME,LONLAT\n');
			for id=1:26
				for t=ts
					if ~isnan(lons_clean(id,t))
						fprintf(fid,'%i,%sT%0.2i%s.000Z,%0.7f %0.7f\n',id,timestamps_str(t,1:10),str2num(timestamps_str(t,12:13))+3,timestamps_str(t,14:19),lons_clean(id,t),lats_clean(id,t));
					end
				end
			end
			fclose(fid);
		end
	end
end	