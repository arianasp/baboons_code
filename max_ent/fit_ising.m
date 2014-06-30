function [ alpha_fit, beta_fit ] = fit_ising( data, weight, converge_thresh, max_iter, init_alphas, init_betas )
%Fits an ising model to the data.
%INPUTS:
%   data: [N x n_samps matrix] of data
%   weight: [number] max value to change parameters by at each step
%   converge_thresh: minimum change in KL divergence to consider the model
%       converged
%   max_iter: [number] maximum number of iterations to do before exiting
%   init_alphas: [N x 1 vector] of alphas to initialize fit
%   init_betas: [N x N matrix] of betas to initialize fit
%OUTPUTS:
%   alphas: [N x 1 vector] of fitted fields
%   betas: [N x N matrix] of fitted couplings

N = size(data,1); %number of individuals
n_samps = size(data,2); %number of samples

%default parameters
if ~exist('init_alphas')
    init_alphas = rand(1,N) / (N + (N^2 - N)/2); %default to Oren's initialization method (random nums between 0 and 1, divide by number of parameters)
end

if ~exist('init_betas')
    init_betas = rand(N,N) / (N + (N^2 - N)/2); %default to Oren's initialization method
    for i = 1:N
        for j = (i+1):N
            init_betas(j,i) = init_betas(i,j);
        end
        init_betas(i,i) = 0;
    end
end
    
if ~exist('max_iter')
    max_iter = 100000;
end

if ~exist('converge_thresh')
    converge_thresh = 10^(-7);
end

if ~exist('weight')
    weight = 0.5;
end

%compute state matrix
[ state_mat ] = generate_state_mat( N );

%get empirical probabilities of each state
[ P_emp ] = empirical_state_probabilities( data );

%get expectations for xi and xi,xj for the data
E_emp = state_mat * P_emp;

tri_mat = triu(ones(N,N),1); %generate upper triangular matrix to use for indexing into couplings matrix

init_coeffs = [init_alphas init_betas(find(tri_mat==1))']; %initialize coefficient vector

coeffs = init_coeffs; %initialize coefficients

%run the fit
previous_KL_div = Inf;
for iter = 1:max_iter
    
    P_model = ising_probabilities( coeffs, state_mat ); %get the model probabilities
    
    E_model = state_mat * P_model; %get the model expectation values
    
    E_diff = E_emp - E_model; %difference between model and data expectation values
    
    tot_E_diff = sum(abs(E_diff)); %total difference in expectations (used as convergence condition)
    
    curr_KL_div = KL_divergence(P_emp,P_model); %compute KL divergence btwn empirical & current model

    dKL = previous_KL_div - curr_KL_div; %get the change in KL divergence
    
    %if converged
    if dKL < converge_thresh
        alpha_fit = coeffs(1:N);
        beta_fit = tri_mat;
        beta_fit(find(tri_mat==1))=coeffs((N+1):end);
        for i = 1:N
            for j = (i+1):N
                beta_fit(j,i) = beta_fit(i,j);
            end
        end
        fprintf('Model converged in %d iterations. \n',iter)
        fprintf('KL = %f \n',curr_KL_div)
        fprintf('dKL = %f \n',dKL)
        
        return 
    %if not converged    
    else
        curr_weight = weight; %initialize current weight to its maximum value ("weight")
        new_KL_div = Inf; %initialize new KL divergence at infinity
        
        while new_KL_div > curr_KL_div
            new_coeffs = coeffs + curr_weight*E_diff'; %get new coefficients using current weight
            new_P_model = ising_probabilities( new_coeffs, state_mat ); %get probabilities from new model
            new_KL_div = KL_divergence(P_emp, new_P_model); %KL divergence btwn empirical and new model
            curr_weight = curr_weight / 2; %current weight gets divided in half, until new model is better than old (according to KL divergence)
        end

        coeffs = coeffs + curr_weight*E_diff'; %update coefficients
        previous_KL_div = curr_KL_div; %update previous KL divergence 
        
    end
end

fprintf('Model failed to converge in %d iterations \n',iter)
fprintf('KL = %f \n',curr_KL_div)
fprintf('dKL = %f \n',dKL)
warning('model failed to converge - returning current coefficients')
alpha_fit = coeffs(1:N);
beta_fit = tri_mat;
beta_fit(find(tri_mat==1))=coeffs((N+1):end);
for i = 1:N
	for j = (i+1):N
        beta_fit(j,i) = beta_fit(i,j);
	end
end




end

