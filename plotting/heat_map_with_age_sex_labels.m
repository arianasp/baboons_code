function [ output_args ] = heat_map_with_age_sex_labels( data_mat, baboon_info, clims, ind_ranks )
%Creates a heat map from a matrix, and also puts colored markers to
%indicate the axe/sex class of each individual
%INPUTS:
%   data_mat: [NxN matrix] of values representing some kind of dyadic
%       interaction between N individuals (e.g. leader/follower relation)
%   baboon_info: [struct] containing age / sex class and other information
%   clims: [2x1 vector] color scale range (defaults to [min_val max_val])
%   ind_ranks: [Nx1 vector] giving the row/col that each individual in the
%       baboon_info struct occupes in the data matrix
%OUTPUTS:
%   a plot

%default clims
if ~exist('clims')
    clims = [nanmin(nanmin(data_mat)) nanmax(nanmax(data_mat))];
end

%number of individuals
N = length(baboon_info);

if size(data_mat) ~= [N N];
    error('data matrix must be NxN, where N is also the length of the baboon_info struct')
end

figure
hold on;
imagesc(data_mat,clims)
colormap(red_white_blue_colormap())
axis square
axis([-3 N+1 -3 N+1])

for i = 1:N
    y = ind_ranks(i);
    
    %get appropriate color for age/sex class
    if strcmp(baboon_info(i).sex,'M')
    	if strcmp(baboon_info(i).age,'A')
        	col = [0 0 255]./255;
        elseif strcmp(baboon_info(i).age,'SA')
        	col = [102 204 255]./255;
        else
        	col = [102 102 102]./255;
        end
    else
    	if strcmp(baboon_info(i).age,'A')
    		col = [255 0 0]./255;
    	elseif strcmp(baboon_info(i).age,'SA')
    		col = [255 204 102]./255;
    	else
    		col = [102 102 102]./255;
    	end
   end
    
    %get appropriate marker shape for sex
    if strcmp(baboon_info(i).sex,'M')
        mark = 'o';
    elseif strcmp(baboon_info(i).sex,'F')
        mark = 'o';
    end
    
    plot(0,y,mark,'MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor',col)
    plot(y,0,mark,'MarkerSize',10,'MarkerFaceColor',col,'MarkerEdgeColor',col)
    h1=text(-1.5,y,baboon_info(i).collar_num,'Color','black','FontWeight','bold','HorizontalAlignment','Center');
    h2=text(y,-1.5,baboon_info(i).collar_num,'Color','black','FontWeight','bold','HorizontalAlignment','Center');
	set(h2,'rotation',90)
end
    
    
    
end

