function [ H ] = indep_entropy( states )
%Compute the entropy of an independent distribution

N = size(states,1);

means = mean(states,2)';
indep_probs = independent_probabilities(means);

H_i = indep_probs .* log2(indep_probs);
H_i(find(indep_probs == 0)) = 0;

H=-sum(H_i);


end

