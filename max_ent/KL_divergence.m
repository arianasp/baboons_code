function [ KL_div ] = KL_divergence( P, Q )
%Compute the KL divergence of two probability distributions (defined as two
%vectors containing probabilities). D_KL(P||Q)
%INPUTS:
%   P: [N x 1 vector] of probabilities of each state in distribution P
%   Q: [N x 1 vector] of probabilities of each state in distribution Q

PlogP = P.*log(P);
PlogP(find(P == 0)) = 0;
PlogQ = P.*log(Q);
KL_div = sum(PlogP - PlogQ);

end

