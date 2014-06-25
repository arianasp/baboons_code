function [ E1, E2 ] = expectations_from_state_probabilities( probs, state_mat )
%INPUTS:
%   probs: [2^N x 1 vector] of state probabilities
%   state_mat: [(N+N^2) x 2^N matrix] where each row represents a possible
%       state (see generate_state_mat for details)
%OUTPUTS:
%   E1: expectation for the value of a given individual
%   E2: expeectation for the value of a pair

N = log2(size(probs,1));




end

