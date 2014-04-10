function [ DS ] = davids_score( dom_mat )
%Computes David's Score from a matrix of dominance interactions between
%pairs of individuals. The DS is defined as in Balasubramaniam et al. 2013
% DS = (W + W2) - (L + L2)
%INPUTS:
%   dom_mat: [NxN matrix] of dominance interactions, where 
%       dom_mat(i,j) = number of times in which i dominates j. 
%OUTPUTS:
%   DS: [N vector] of the David's score for each individual

N = size(dom_mat,1); %number of individuals

P = zeros(N,N); %proportion of wins for each pair
L = zeros(N,N); %proportion of losses for each pair
for i = 1:N
    for j = 1:N
        if not(dom_mat(i,j) + dom_mat(j,i) == 0 || i==j)
            P(i,j) = (dom_mat(i,j)) / (dom_mat(i,j) + dom_mat(j,i));
            L(i,j) = (dom_mat(j,i)) / (dom_mat(i,j) + dom_mat(j,i));
        end
    end
end


w = sum(P,2);
l = sum(L,2);

w_weighted = zeros(N,1);
l_weighted = zeros(N,1);
for i = 1:N
    pw = 0;
    lw = 0;
    for j = 1:N
        if not(i==j)
            pw = pw + w(j)*P(i,j);
            lw = lw + l(j)*L(i,j);
        end
    end
    w_weighted(i) = pw;
    l_weighted(i) = lw;
end

P
L

[w w_weighted l l_weighted]

DS = w + w_weighted - l - l_weighted;


 


end

