function [ state_probs ] = independent_probabilities( indep_probs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


N = length(indep_probs);

state_probs = zeros(1,2^N);

one_minus_probs = ones(1,N) - indep_probs;

for idx = 1:(2^N)
    xi = dec2bin(idx-1,N) - '0';
    
    up = find(xi==1);
    down = find(xi==0);
    
    state_probs(idx) = prod(indep_probs(up)) * prod(one_minus_probs(down));

    
end



end

