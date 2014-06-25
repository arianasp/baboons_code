function [ data ] = generate_ising_data( alpha, beta, n_samps )
%Generates data from an ising model with parameters alpha (Nx1 vector of 
%fields) and beta (NxN matrix of couplings).
%INPUTS:
%   alpha: [N x 1 vector] of biases / fields 
%   beta: [N x N matrix] of couplings
%   n_samps: [number] of samples to generate using the ising model
%OUTPUTS:
%   data =  [N x n_samps] matrix of generated data

N = length(alpha); %number of individuals

state_mat = generate_state_mat(N); %generate state matrix for fast computation

tri_mat = triu(ones(N,N),1); %generate upper triangular matrix to use for indexing into couplings matrix

coeffs = [alpha beta(find(tri_mat==1))'];

%compute the probability of each state
[ probs ] = ising_probabilities( coeffs, state_mat );

cumprobs = cumsum(probs); %cumulative probability distribution (for drawing from)

data = zeros(N,n_samps);
for i = 1:n_samps
    p = rand(1); %pick a random number
    
    idx = find(p < cumprobs,1,'first') -1; %find the associated state
    
    curr_dat = dec2bin(idx,N) - '0'; %convert to binary
    
    data(:,i) = curr_dat; %store data
end


end

