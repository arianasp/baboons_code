function [ output_args ] = rank_visualization( ranks, baboon_info, ori, context1, context2)
%Make a visualization of the ranks of each baboon in different context.
%INPUTS:
%   ranks: [NxM matrix] where N is the number of individuals and M is the
%       number of contexts in which they are ranked.
%   baboon_info: [length N struct] of baboon information
%   ori: whether to orient things vertical ('vert' or 'vertical') or
%       horizontally ('horiz' or 'horizontal')
%   context1: [string] name of context 1
%   context2: [string] name of context 2

N = size(ranks,1);
M = size(ranks,2);

if length(baboon_info) ~= N
    error('number of rows in "ranks" must equal length of "baboon_info"')
end

colors = [255 0 0; 255 204 102; 0 0 255; 102 204 255; 102 102 102] ./ 255;

h=figure;
set(h,'Position',[100 100 1200 600])
set(h,'PaperPositionMode','auto')
hold on;
for i = 1:N
    
    if strcmp(baboon_info(i).sex,'F')
        if strcmp(baboon_info(i).age,'A')
            color = colors(1,:);
        elseif strcmp(baboon_info(i).age,'SA')
            color = colors(2,:);
        else
            color = colors(5,:);
        end
    else
        if strcmp(baboon_info(i).age,'A')
            color = colors(3,:);
        elseif strcmp(baboon_info(i).age,'SA')
            color = colors(4,:);
        else
            color = colors(5,:);
        end
    end

    
    if strcmp(ori,'vert') || strcmp(ori,'vertical')
        plot(1:M,ranks(i,:),'-','LineWidth',3,'Color',color)
        plot(1:M,ranks(i,:),'.','MarkerSize',30,'Color',color)
        text(1-.1,ranks(i,1),baboon_info(i).collar_num,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
        text(2+.1,ranks(i,2),baboon_info(i).collar_num,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
    elseif strcmp(ori,'horiz') || strcmp(ori,'horizontal')
        plot(ranks(i,:),1:M,'-','LineWidth',3,'Color',color)
        plot(ranks(i,:),1:M,'.','MarkerSize',30,'Color',color)
        text(ranks(i,1),1-.1,baboon_info(i).collar_num,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
        text(ranks(i,2),2+.1,baboon_info(i).collar_num,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
    end
end

if strcmp(ori,'horiz') || strcmp(ori,'horizontal')
    axis([0 N+1 0 M+1])
else
    axis([0 M+1 0 N+1])
end

if exist('context1') && exist('context2')
    if strcmp(ori,'horiz') || strcmp(ori,'horizontal')
        text(length(ranks)/2,.7,context1,'FontSize',24,'HorizontalAlignment','center','FontWeight','bold')
        text(length(ranks)/2,2.3,context2,'FontSize',24,'HorizontalAlignment','center','FontWeight','bold')
    end
end




end

