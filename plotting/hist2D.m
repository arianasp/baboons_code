function [ Z ] = hist2D( data, nbins, xrange, yrange )
%Makes a 2d histogram (heat map) of the data in [x y] 
%INPUTS:
%   data: [N x 2 (or 2 x N) matrix] of N data points in 2 dimensions
%   nbins: [2 x 1 vector] of the number of bins to use in each dimension
%   xrange: [2 x 1 vector] of the min and max x values
%   yrange: [2 x 1 vector] of the min and max y values
%OUTPUTS:
%   a 2d histogram (heat map)

if size(data,1) ~= 2 && size(data,2) ~= 2
    error('data matrix must be Nx2')
end

if size(data,1) == 2
    data = data';
end

xi = transpose(linspace(xrange(1), xrange(2),nbins(1)));
yi = transpose(linspace(yrange(1), yrange(2),nbins(2)));

xr = interp1(xi,1:numel(xi),data(:,1),'nearest');
yr = interp1(yi,1:numel(yi),data(:,2),'nearest');

Z = accumarray([xr,yr],1,nbins);

end

