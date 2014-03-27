function [ ] = ind_hists_by_age_sex_class( data, baboon_info, day_start_idxs, day_range, bins )
%Plots a (normalized) histogram of the values in data, with all individuals
%plotted on the same axes (with the same bins). 
%INPUTS:
%   data: [NxT matrix] of data values for each individuals as a function of
%       time
%   baboon_info: [struct length N] containing age/sex class info among 
%       other things
%   day_start_idxs: [n_days x 1 vector] containing the indexes to the start
%       of each day
%   day_range: [vector] indicating which days to use (e.g. 1:14)
%   bins: [vector] indicating the bins to use for the histograms
%OUTPUTS:
%   a histogram with lines for each individual, colored according to
%   age/sex class. black = adult, red = subadult, green = juvenile. solid =
%   male, dotted = female

%get number of individuals at time steps
N = size(data,1);
T = size(data,2);

%add an ending time to day_start_idxs
day_start_idxs = [day_start_idxs T+1];

%get the appropriate data (over all days in day_range)
hist_data = [];
for d = day_range
    hist_data = [hist_data data(:,day_start_idxs(d):(day_start_idxs(d+1)-1))];
end

%make the figure
figure
hold on;
for i = 1:N
    histo = hist(hist_data(i,:),bins);
    histo = histo ./ sum(histo);
    if strcmp(baboon_info(i).age,'A')
        col = 'black';
    elseif strcmp(baboon_info(i).age,'SA')
        col = 'red';
    elseif strcmp(baboon_info(i).age,'J')
        col = 'green';
    else
        error('unknown age class')
    end
    
    if strcmp(baboon_info(i).sex,'M')
        line = '-';
    elseif strcmp(baboon_info(i).sex,'F')
        line = '--';
    else
        error('unknown sex class')
    end
    
    plot(bins,histo,'LineStyle',line,'Color',col,'LineWidth',2)
end


end

