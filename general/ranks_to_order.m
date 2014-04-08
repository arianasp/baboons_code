function [ order ] = ranks_to_order( ranks )
%Takes a vector giving the ranks of individuals and outputs a vector of the
%order. In other words, "order" gives the order of the
%individual (e.g. 5 first, 1 second, 20 third = [5,1,20...]) and "ranks" gives the
%rankings of each individual in that order (e.g. 1 was ranked 2nd, 2 was
%ranked 15th, 3 was ranked 7th = [2 15 7 ...]).

N = length(ranks);

order = transpose(1:N);

order(ranks) = order;

end

