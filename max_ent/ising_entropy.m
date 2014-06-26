function [ H ] = ising_entropy( alpha, beta )
%Compute the entropy of the distribution specified by an ising model with
%parameters alpha and beta.
%INPUTS:
%   alpha: [N x 1 vector] of alpha values (fields / biases)
%   beta: [N x N matrix] of beta values (couplings between individuals)
%OUTPUTS:
%   E: entropy of the distribution

N = length(alpha);

state_mat = generate_state_mat(N);

tri_mat = triu(ones(N,N),1); %generate upper triangular matrix to use for indexing into couplings matrix

coeffs = [alpha beta(find(tri_mat==1))']; %coefficient vector

[ probs ] = ising_probabilities( coeffs, state_mat );

H_i = probs .* log2(probs);

H_i(find(probs == 0)) = 0;

H = -sum(H_i);


end

