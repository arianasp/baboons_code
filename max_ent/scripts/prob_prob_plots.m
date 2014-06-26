%Prob/prob plots

%PARAMS
filename = '/Users/arianasp/Desktop/Baboons/output/max_ent/baboon_ising_fit2.mat';
savename = '/Users/arianasp/Desktop/Baboons/output/max_ent/plots/N15_14days_prob_prob_plots.png';

%load data
load(filename)
N = size(baboon_states,1);
nsamp = size(baboon_states,2);

%COMPUE PROB DISTRIBUTIONS

%empirical
[ P_data ] = empirical_state_probabilities( baboon_states );

%independent
indep_probs = mean(baboon_states,2);
[ P_indep ] = independent_probabilities( indep_probs' );

%ising
[ state_mat ] = generate_state_mat( N );
tri_mat = triu(ones(N,N),1); %generate upper triangular matrix to use for indexing into couplings matrix
coeffs = [alpha_fit beta_fit(find(tri_mat==1))']; %initialize coefficient vector
[ P_ising ] = ising_probabilities( coeffs, state_mat );

%MAKE PLOTS
figure

%independent vs data
subplot(1,2,1)
hold on;
plot(P_data,P_indep,'.')
axis([1/nsamp 0.06 1/nsamp 0.06])
set(gca,'XScale','log','YScale','log')
loglog([1/nsamp 0.06],[1/nsamp 0.06],'-k')
axis square
xlabel('P_{data}')
ylabel('P_{indep}')

subplot(1,2,2)
hold on;
plot(P_data,P_ising,'.')
axis([1/nsamp 0.06 1/nsamp 0.06])
set(gca,'XScale','log','YScale','log')
loglog([1/nsamp 0.06],[1/nsamp 0.06],'-k')
axis square
xlabel('P_{data}')
ylabel('P_{ising}')

print('-dpng',savename)


