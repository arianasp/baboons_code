function [ output_args ] = dyadic_interaction_movie( xs, ys, interactions, idx, save_imgs, outdir )
%Makes a movie of a dyadic interaction between a pair of individuals.
%INPUTS:
%   xs: [NxT matrix] x positions for all individuals
%   ys: [NxT matrix] y positions for all individuals
%   interactions: [struct] containing interactions between a pair of
%       individuals (output from dyadic_interactions.m)
%   idx: [number] which interaction to make a movie of
%   save_imgs: whether to save the images as jpegs
%   outdir: the directory in which to save the images
%OUTPUTS:
%   No outputs, but pops up images of where the individuals are at each
%       time point during the event, and saves them if requested. The
%       "leader" is colored blue and the "follower" is colored red, all
%       other individuals are colored black. A line is drawn between the
%       leader and the follower for visualization. The line is black during
%       the time period t1 to t2, then turns red for the time period t2 to
%       t3.

if ~exist('save_imgs')
    save_imgs = 0;
end

if ~exist('outdir')
    outdir = '/Users/arianasp/Desktop/Baboons/movies/slinky_interactions';
end

%get leader and follower
lead = interactions(idx).leader;
foll = interactions(idx).follower;

%get times
times = interactions(idx).time_idxs;

t_range = times(1):times(3);
%get size of area needed to plot
minx = nanmin(nanmin(xs(:,t_range)))-20;
miny = nanmin(nanmin(ys(:,t_range)))-20;
maxx = nanmax(nanmax(xs(:,t_range)))+20;
maxy = nanmax(nanmax(ys(:,t_range)))+20;


%set up save location, if required
if save_imgs
    dir = [outdir '/' interactions(idx).type '_strength' num2str(round(interactions(idx).strength_mult)) '_disp' num2str(round(interactions(idx).disp_mult)) '_lead' num2str(interactions(idx).leader) '_foll' num2str(interactions(idx).follower) '_t' num2str(times(1)) '-' num2str(times(3))];
    mkdir(dir)
end
    

figure
i = 1;
for t = t_range
    hold on;
    x = xs(:,t);
    y = ys(:,t);
    
    if t <= times(2)
        plot([x(lead) x(foll)],[y(lead) y(foll)],'black')
    else
        plot([x(lead) x(foll)],[y(lead) y(foll)],'red')
    end
    
    plot(x,y,'o','MarkerSize',5,'MarkerEdgeColor','black','MarkerFaceColor','black')
    plot(x(lead),y(lead),'o','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','blue')
    plot(x(foll),y(foll),'o','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red')
    axis([minx maxx miny maxy])
    axis equal
    title_str = [interactions(idx).type ', strength = ' num2str(interactions(idx).strength_mult) ', disparity = ' num2str(interactions(idx).disp_mult)];
    title(title_str)
    xlabel(['t = ' num2str(t)])
    if save_imgs
        print('-dpng',[dir '/' num2str(i) '.png'])
    else
        pause(0.05)
    end
    i=i+1;
    
    clf
end


end

