function [ color_vector ] = red_white_blue_colormap( )
%Outputs a colormap that goes from blue --> white and then white --> blue
%in 256 steps
%INPUTS:
%   none
%OUTPUTS:
%   [256x3 vector] of color triplets

color_vector = ones(256,3);

color_vector(1:128,1) = 0:(1/127):1;
color_vector(1:128,2) = 0:(1/127):1;
color_vector(129:256,2) = 1:(-1/127):0;
color_vector(129:256,3) = 1:(-1/127):0;


end

