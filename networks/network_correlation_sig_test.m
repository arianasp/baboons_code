function [ stats ] = network_correlation_sig_test( A1, A2, corr_type, net_type, n_randomizations )
%Computes the network correlation and performs a randomization test to test
%the significance of the result.
%INPUTS:
%   A1: [N x N matrix] of network connections 
%   A2: [N x N matrix] of network connections
%   corr_type: [string] type of correlation to use. Right now, supported
%   types are: 'Pearson'
%   net_tye: [string] whether the network is 'directed' (in which case the 
%       correlation is run over all entries) or 'undirected'/'symmetric', 
%       or 'antisymmetric' (in which case the correlation is run only over
%       the upper triangular entries)
%OUTPUTS:
%   stats: [struct] with fields
%       .R: value of the correlation
%       .p: p-value

%get correlation between the two networks
[ R_data ] = network_correlation( A1, A2, corr_type, net_type );

R_null = zeros(n_randomizations,1);
for r = 1:n_randomizations
    [ A1_perm ] = permute_edges( A1, net_type ); %permute one of the matrices
    R_null(r) = network_correlation( A1_perm, A2, corr_type, net_type );
end

p = sum(R_null < R_data) / n_randomizations;

stats.R = R_data;
stats.p = p;


end

