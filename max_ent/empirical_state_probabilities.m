function [ probs ] = empirical_state_probabilities( states )
%Computes the probability of each state from a matrix of data.
%INPUTS:
%   states: [N x n_samp matrix] of data - the states of N individuals from
%       n_samp samples. Should contain binary states (0s and 1s)
%OUTPUTS:
%   probs: [2^N x 1 vector] of the empirical probability of each state
%NOTE:
%   This code only considers bianry states.

N = size(states,1); %number of individuals
n_samp = size(states,2); %number of samples

if n_samp < N
    error('number of samples is less than number of individuals - you may need to transpose your input matrix')
end

probs = zeros(2^N,1); %probability of each state occurring in the data

for i = 1:n_samp
    currdat = states(:,i)'; %get the current sample
    idx = sum(currdat.*2.^(numel(currdat)-1:-1:0)) + 1; %convert to decimal
    probs(idx) = probs(idx) + 1; %add to probability vector
end

probs = probs / sum(probs); %normalize probabilities



end

