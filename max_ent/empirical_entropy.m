function [ H ] = empirical_entropy( states, epsilon )
%Computes the entropy from an empirical list of observed states
%INPUTS:
%   states: [N x n_samp vector] of states
%OUTPUTS;
%   H: entropy of the distribution


N = size(states,1);
n_samp = size(states,2);

[ probs ] = empirical_state_probabilities(states);

%add-one rule (actually, add-epsilon rule), this is to get rid of zeros
%probs = probs + (epsilon / n_samp);
%probs = probs / sum(probs); %renormalize

H_i = probs .* log2(probs);

H_i(probs == 0) = 0;

H = -sum(H_i);

end

