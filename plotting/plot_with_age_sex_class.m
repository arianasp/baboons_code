function [ output_args ] = plot_with_age_sex_class( data_x, data_y, baboon_info, x_label, y_label, plot_title )
%Plots the data in the vectors data_x and data_y with colors and shapes
%representing the age sex class of individual.
%INPUTS:
%   data_x : [Nx1 vector] x data to plot 
%   data_y : [Nx1 vector] y data to plot [Nx1 vector]
%   baboon_info: [struct] containing age/sex class information (among other
%   things)
%OUTPUTS:
%   a plot of data_y vs data_x, with the following markers:
%       males: squares
%       females: circles
%       adults: black
%       subadults: red
%       juveniles: green

N = length(baboon_info);

if length(data_x)~= N || length(data_y) ~= N
    error('data vectors must be the same length as baboon_info struct')
end

adults = nan(N,1); subadults = nan(N,1); juveniles = nan(N,1);
males = nan(N,1); females = nan(N,1);

for i = 1:N
    
    if strcmp(baboon_info(i).age,'A')
        adults(i)=1;
    elseif strcmp(baboon_info(i).age,'SA')
        subadults(i)=1;
    elseif strcmp(baboon_info(i).age,'J')
        juveniles(i) = 1;
    else
        error('unknown age class')
    end
    
    if strcmp(baboon_info(i).sex,'M')
        males(i)=1;
    elseif strcmp(baboon_info(i).sex,'F')
        females(i)=1;
    else
        error('unknown sex class')
    end
end

figure
hold on;

set(0,'DefaultAxesFontSize',16)

plot(data_x.*males.*adults,data_y.*males.*adults,'o','MarkerSize',10,'MarkerFaceColor',[0 0 255]./255,'MarkerEdgeColor',[0 0 255]./255)
plot(data_x.*males.*subadults,data_y.*males.*subadults,'o','MarkerSize',10,'MarkerFaceColor',[102 204 255]./255,'MarkerEdgeColor',[102 204 255]./255)
plot(data_x.*males.*juveniles,data_y.*males.*juveniles,'o','MarkerSize',10,'MarkerFaceColor',[102 102 102]./255,'MarkerEdgeColor',[102 102 102]./255)
plot(data_x.*females.*adults,data_y.*females.*adults,'o','MarkerSize',10,'MarkerFaceColor',[255 0 0]./255,'MarkerEdgeColor',[255 0 0]./255)
plot(data_x.*females.*subadults,data_y.*females.*subadults,'o','MarkerSize',10,'MarkerFaceColor',[255 204 102]./255,'MarkerEdgeColor',[255 204 102]./255)
plot(data_x.*females.*juveniles,data_y.*females.*juveniles,'o','MarkerSize',10,'MarkerFaceColor',[102 102 102]./255,'MarkerEdgeColor',[102 102 102]./255)

xlabel(x_label)
ylabel(y_label)
title(plot_title)


end

