function [ H_ests ] = entropy_vs_sample_size( states, data_fracs, reps )
%Estimates the entropy of a distribution as a function of the number of
%samples taken, using jackknifing. 
%INPUTS:
%   states: [N x nsamp matrix] of states
%   data_fracs: [vector] giving the inverse fraction of the data to use in 
%       making each estimate (defaults to 1:10)
%   reps: [number] giving the number of times to repeat the subsampling and
%       compute entropy (defaults to 1)
%OUTPUTS:
%   H_ests: [n_fracs x reps matrix] giving the estimated entropy for each 
%       fraction, for each replicat of the jackknifing procedure

if ~exist('data_fracs')
    data_fracs = 1:10;
end

if ~exist('reps')
    reps = 1;
end

N = size(states,1);
nsamp = size(states,2);

H_ests = zeros(length(data_fracs),reps);
for r = 1:reps
    for f = 1:length(data_fracs);
        f
        frac = data_fracs(f);
        idxs = randsample(nsamp,round(nsamp/frac));
        states_sub = states(:,idxs);
        [ H ] = empirical_entropy( states_sub);
        H_ests(f,r)=H;
    end
end



end

