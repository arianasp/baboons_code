%Test ising model fitting by generating data from a given set of parameters
%and then fitting the model, then checking whether the parameters agree

%INPUT PARAMETERS HERE
N = 14; %number of individuals
n_samps = 100000; %number of samples to take
weight = 1;
converge_thresh = 0.001;
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
disp('generating data')
tic
[ data ] = generate_ising_data( alpha, beta, n_samps );
toc

%Fit an Ising model to the generated data (use default parameters in the
%fit)
disp('fitting ising model')
tic
[ alpha_fit, beta_fit ] = fit_ising( data, weight, converge_thresh, max_iter);
disp('model fit!')
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




