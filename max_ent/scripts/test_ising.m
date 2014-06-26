%Test ising model fitting by generating data from a given set of parameters
%and then fitting the model, then checking whether the parameters agree

%INPUT PARAMETERS HERE
N = 10; %number of individuals
n_samps = 100000; %number of samples to take
weight = 0.1;
converge_thresh = 10^(-7);
max_iter = 100000;

%Set some random parameters (between -1 and 1)
alpha = (rand(1,N)-.5)*2;
beta = (rand(N,N)-.5)*2;
for i = 1:N
    for j = (i+1):N
        beta(j,i) = beta(i,j);
    end
    beta(i,i) = 0;
end

%Generate data with these parameters
disp('Generating data...')
tic
[ data ] = generate_ising_data( alpha, beta, n_samps );
toc

%Fit an Ising model to the generated data (use default parameters in the
%fit)
disp('Fitting ising model...')
tic
[ alpha_fit, beta_fit ] = fit_ising( data, weight, converge_thresh, max_iter);
toc

%Plot results
figure

subplot(2,2,1)
imagesc(beta,[-1 1])
title('beta values')
axis square
colormap(red_white_blue_colormap())
colorbar

subplot(2,2,2)
imagesc(beta_fit,[-1 1])
title('fitted beta values')
axis square
colormap(red_white_blue_colormap())
colorbar

subplot(2,2,3)
hold on
plot([-1 1],[-1 1],'-k')
plot(alpha,alpha_fit,'.')
xlabel('alpha')
ylabel('fitted alpha')
axis([-1 1 -1 1])

subplot(2,2,4)
imagesc(beta_fit - beta,[-1 1])
title('difference between beta and fitted beta')
axis square
colormap(red_white_blue_colormap())
colorbar

%Compute entropy
H_orig = ising_entropy(alpha,beta); %compute entropy of model distribution
H_sim = empirical_entropy(data,0); %entropy of simulated data
H_fit = ising_entropy(alpha_fit, beta_fit); %entropy of the fitted model

fprintf('Entropy of the original ising model: %f \n',H_orig);
fprintf('Entropy of the simulated data distribution: %f \n',H_sim);
fprintf('Entropy of the fitted ising model: %f \n',H_fit);




