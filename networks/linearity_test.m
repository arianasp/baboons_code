function [ p, ntriads_data, ntriads_null ] = linearity_test( adj_mat, n_rand )
%Tests whether a network (given by the directed adjacency matrix adj_mat)
%is more 'linear' (aka hierarchical) than random. Does this using the
%double randomization method of de Vries 1995, which considers only triads 
%and not higher-order cycles. adj_mat(i,j) should be 1 if i leads j, 0 if j
%leads i, and  NaN if no directional relationship can be determined
%INPUTS:
%   adj_mat: [NxN matrix] adjacency matrix
%OUTPUT:
%   p: [number] p-value for testing the hypothesis that adj_mat contains
%       fewer cyclical triads (i.e. is more linear/hierarchical) than 
%       expected based on a randomized null model
%   ntriads_data: [n_rand x 1 vector] of the number of cyclical triads 
%       present in adj_mat
%       (the reason this is a list of values instead of a single number is
%       that when no relationship can be determined between two nodes, one
%       nodes is randomly assigned to point to the other by the algorithm
%       during each iteration
%   ntriads_null: [n_rand x 1 vector] list of the number of cyclical triads
%       present in each randomized null model


%remove rows and columns that are only NaNs
i = 1;
while i <= size(adj_mat,1)
    row = adj_mat(i,:);
    col = adj_mat(:,i);
    if sum(not(isnan(row)))==0
        adj_mat = [adj_mat(1:(i-1),:) ; adj_mat((i+1):end,:)];
        adj_mat = [adj_mat(:,1:(i-1))   adj_mat(:,(i+1):end)];
        i = i - 1;
    end
    i = i + 1;
end

d_sig_count = 0;
ntriads_data = zeros(n_rand,1);
ntriads_null = zeros(n_rand,1);
for i = 1:n_rand
    
    %randomize missing links
    adj_mat_real = randomize_directional_links(adj_mat,'missing');
    d_real = compute_n_triads(adj_mat_real);
    ntraids_data(i) = d_real;
    
    %randomize existing links
    adj_mat_fake = randomize_directional_links(adj_mat,'existing');
    adj_mat_fake(find(isnan(adj_mat))) = adj_mat_real(find(isnan(adj_mat)));
    d_fake = compute_n_triads(adj_mat_fake);
    ntriads_null(i) = d_fake;
    
    if d_fake <= d_real
        d_sig_count = d_sig_count + 1;
    end
end

p = d_sig_count / n_rand;




end

