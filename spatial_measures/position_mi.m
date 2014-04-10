function [ MI ] = position_mi( xs, ys, i, j, time_range, xmin, xmax, ymin, ymax, nbins_x, nbins_y )
%Computes the mutual information between the positions of two
%individuals during the time range time_range;
%Uses the function mi from MIToolbox for the main mutual information
%computation.

xi = xs(i,time_range);
xj = xs(j,time_range);
yi = ys(i,time_range);
yj = ys(j,time_range);

[ subsi ] = spatial_bin( xi,yi,xmin,xmax,ymin,ymax,nbins_x,nbins_y  );
[ subsj ] = spatial_bin( xj,yj,xmin,xmax,ymin,ymax,nbins_x,nbins_y  );

narows = find(sum(isnan(xi),1)>0);
narows = intersect(narows,find(sum(isnan(xj),1)>0)); 
narows = intersect(narows,find(sum(isnan(yi),1)>0)); 
narows = intersect(narows,find(sum(isnan(yj),1)>0)); 

narows

size(subsi)
size(subsj)

subsi(narows,:) = [];
subsj(narows,:) = [];

MI = mi(subsi,subsj);



end

