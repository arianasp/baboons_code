function [ ranks ] = order_to_ranks( order )
%Takes a vector giving the order of individuals in some ranking and outputs
%a vector of their ranks. In other words, "order" gives the order of the
%individual (e.g. 5 first, 1 second, 20 third = [5,1,20...]) and "ranks" gives the
%rankings of each individual in that order (e.g. 1 was ranked 2nd, 2 was
%ranked 15th, 3 was ranked 7th = [2 15 7 ...]).

N = length(order);

ranks = order;

ranks(ranks) = 1:N;

end

