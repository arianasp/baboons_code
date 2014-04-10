function [ subs ] = spatial_bin( x,y,xmin,xmax,ymin,ymax,nbins_x,nbins_y  )
%Bin spatial data to produce discrete values from continuous ones. Takes in
%two vectors, x and y, of spatial data points (x,y) and puts them into
%rectangular bins. x values range from xmin to xmax, and there are nbins_x
%bins in the x direction (and correspondingly for y bins). Values (x,y) 
%outside of the range xmin < x < xmax and ymin < y < ymax are associated
%with the nearest bin.
%This function preserves NaNs.
%INPUTS:
%   x,


xi = transpose(linspace(xmin, xmax,nbins_x));
yi = transpose(linspace(ymin, ymax,nbins_y));

xr = interp1(xi,1:numel(xi),x,'nearest');
yr = interp1(yi,1:numel(yi),y,'nearest');

size(xr)

if size(xr,2) == 1
    subs = transpose([xr; yr]);
elseif size(xr,1) == 1
    subs = [xr', yr'];
end


end

