%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This file demonstrates the use of the program MultiSpec(), which
%   implements the method discussed in 
%   "Adaptive Bayesian Power Spectrum Analysis of Multivariate Nonstationary Time Series."
%
%   Content:
%   (1) Bivariate slowly-varying time series  (section 4.2) 
%       (1a) Simulate data for bivariate slowly-varying time series
%       (1b) Set the options and hyperparameters for the sampler
%       (1c) Run the Bayesian estimation process
%       (1d) Obtain and plot the number and location of partitions
%       (1e) Obtain and plot the surface of spectra and coherence
%       (1f) Obtain and plot the credible intervals for surfaces
%   (2) Trivariate abrupt-changing time series (section 4.1)
%       (1a) Simulate data for trivariate abrupt-changing time series
%       (1b) Set the options and hyperparameters for the sampler
%       (1c) Run the Bayesian estimation process
%       (1d) Obtain and plot the number and location of partitions
%       (1e) Obtain and plot the surface of spectra and coherence
%       (1f) Obtain and plot the credible intervals for surfaces

%% (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1a) Simulate slowly-varying example discussed in Section 4.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(26132);
sig =[1 0.2; 0.2 1];
cont=1.112;
zt = simMAslow(1024, cont, sig, 0)';
dim=size(zt); nobs=dim(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1b) Set the options and hyperparameters for the sampler
params = struct('nloop',10000, 'nwarmup',2000, ...
                   'nexp_max',10, 'tmin',60, 'prob_mml',0.8 , 'nbasis',10, ...
                   'tau_up_limit',10^4, 'sigmasqalpha',10^5, 'init',3,...
                   'nfreq',50, 'verb',1, 'convdiag',0, 'ee', 0.1);   
                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1c) Run the Bayesian estimation process 
rng(20170531);
[spect_matrices, freq_hat, fit, fit_diag] = MultiSpect(zt,params);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1d) Obtain and plot the number and location of partitions
[posterior_probability_slow] = MultiSpect_partition(zt, fit, params); 
number_of_partitions = (1:params.nexp_max)';
T=table(number_of_partitions, posterior_probability_slow);
T 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1e) Obtain and plot the estimated spectra and square coherence surfaces
[spect_slow, coh_slow] = MultiSpect_surface(zt, spect_matrices, fit, params); 

figure
subplot(1,3,1); meshc( 1:nobs, freq_hat,real(spect_slow{1})); set(gca,'YDir','reverse');  title('$ f_{11}$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 8])
subplot(1,3,2); meshc( 1:nobs, freq_hat,real(spect_slow{2})); set(gca,'YDir','reverse');  title('$ f_{22}$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 8])
subplot(1,3,3); meshc( 1:nobs, freq_hat,real(coh_slow{1})); set(gca,'YDir','reverse');  title('$ \rho_{21}^2$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1f) Obtain and plot the surface of spectra and coherence with credible intervals 
alphalevel = 0.05; %desired alpha level
[spect_surf_slow, coh_surf_slow] = ...
    MultiSpect_interval(zt, spect_matrices, fit, params,alphalevel); 

figure
subplot(1,3,1);
maxv = real(max(max([spect_surf_slow{1}(:,:,1), spect_surf_slow{1}(:,:,2), spect_surf_slow{1}(:,:,3)])));
minv = real(min(min([spect_surf_slow{1}(:,:,1), spect_surf_slow{1}(:,:,2), spect_surf_slow{1}(:,:,3)])));
mesh(1:nobs, freq_hat, real(spect_surf_slow{1}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([minv, maxv]); title('$f_{11}$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(spect_surf_slow{1}(:,:,2))); 
mesh(1:nobs, freq_hat, real(spect_surf_slow{1}(:,:,3)));
hold off
subplot(1,3,2);
maxv = real(max(max([spect_surf_slow{2}(:,:,1), spect_surf_slow{2}(:,:,2), spect_surf_slow{2}(:,:,3)])));
minv = real(min(min([spect_surf_slow{2}(:,:,1), spect_surf_slow{2}(:,:,2), spect_surf_slow{2}(:,:,3)])));
mesh(1:nobs, freq_hat, real(spect_surf_slow{2}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([minv, maxv]); title('$f_{22}$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(spect_surf_slow{2}(:,:,2))); 
mesh(1:nobs, freq_hat, real(spect_surf_slow{2}(:,:,3)));
hold off
subplot(1,3,3);
mesh(1:nobs, freq_hat, real(coh_surf_slow{1}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([0, 1]); title('$ \rho_{21}^2$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(coh_surf_slow{1}(:,:,2))); 
mesh(1:nobs, freq_hat, real(coh_surf_slow{1}(:,:,3)));
hold off



%% (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2a) Simulate abrupt-change example discussed in Section 4.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ma11 = [0.6, 0 0; 0.2, -0.5, 0; 0.1 0.3 0.4];
ma12 = [0.3, 0 0; 0,  0.3 0 ; 0 0 0];
ma21 = [0.6, 0 0; 0.2, 0.5 0 ; -0.1 -0.3 0.4];
ma22 = [0.3, 0 0; 0, 0.3 0 ; 0 0 0];
sig  = [1 0.6 0.6; 0.6 1 0.6; 0.6 0.6 1];
%rng(86691);
%rng(20679);
rng(26132);
%rng(19880901)
%rng(19880106)
%rng(20170920)
%rng(20170921)
%rng(20170922)
%rng(20170831)
%rng(20170901)
zt=simMAtriv(600,300,ma11,ma12,ma21,ma22,sig,1000)';
dim=size(zt); nobs=dim(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1b) Set the options and hyperparameters for the sampler
params = struct('nloop',12000, 'nwarmup',4000, ...
                   'nexp_max',6, 'tmin',60, 'prob_mml',0.8 , 'nbasis',10, ...
                   'tau_up_limit',10^4, 'sigmasqalpha',10^5, 'init',3,...
                   'nfreq',50, 'verb',1, 'convdiag',0, 'ee', 0.1);  
               
%%% Note: he code can be run using default parameters 
%%% [spect_matrices, freq_hat, fit, fit_diag] = MultiSpect(zt)
                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2c) Run the Bayesian estimation process 
rng(20170531);
[spect_matrices, freq_hat, fit, fit_diag] = MultiSpect(zt,params);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (1d) Obtain and plot the number and location of partitions
[posterior_probability] = MultiSpect_partition(zt, fit, params); 
number_of_partitions = (1:params.nexp_max)';
T=table(number_of_partitions, posterior_probability);
T 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2e) Obtain and plot the estimated spectra and square coherence surfaces
[spect, coh] = MultiSpect_surface(zt, spect_matrices, fit, params); 

figure
subplot(1,3,1); meshc( 1:nobs, freq_hat,real(spect{1})); set(gca,'YDir','reverse');  title('$ f_{11}$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 6])
subplot(1,3,2); meshc( 1:nobs, freq_hat,real(spect{2})); set(gca,'YDir','reverse');  title('$ f_{22}$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 6])
subplot(1,3,3); meshc( 1:nobs, freq_hat,real(spect{3})); set(gca,'YDir','reverse');  title('$ f_{33}$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 6])
figure
subplot(1,3,1); meshc( 1:nobs, freq_hat,real(coh{1})); set(gca,'YDir','reverse');  title('$ \rho_{21}^2$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 1])
subplot(1,3,2); meshc( 1:nobs, freq_hat,real(coh{2})); set(gca,'YDir','reverse');  title('$ \rho_{31}^2$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 1])
subplot(1,3,3); meshc( 1:nobs, freq_hat,real(coh{3})); set(gca,'YDir','reverse');  title('$ \rho_{32}^2$','Interpreter','LaTex'); xlabel('time'); ylabel('freq'); zlim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (2f) Obtain and plot the surface of spectra and coherence with credible intervals 
alphalevel = 0.05; %desired alpha level
[spect_surf, coh_surf] = ...
    MultiSpect_interval(zt, spect_matrices, fit, params,alphalevel); 

figure
subplot(1,3,1);
maxv = real(max(max([spect_surf{1}(:,:,1), spect_surf{1}(:,:,2), spect_surf{1}(:,:,3)])));
minv = real(min(min([spect_surf{1}(:,:,1), spect_surf{1}(:,:,2), spect_surf{1}(:,:,3)])));
mesh(1:nobs, freq_hat, real(spect_surf{1}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([minv, maxv]); title('$f_{11}$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(spect_surf{1}(:,:,2))); 
mesh(1:nobs, freq_hat, real(spect_surf{1}(:,:,3)));
hold off
subplot(1,3,2);
maxv = real(max(max([spect_surf{2}(:,:,1), spect_surf{2}(:,:,2), spect_surf{2}(:,:,3)])));
minv = real(min(min([spect_surf{2}(:,:,1), spect_surf{2}(:,:,2), spect_surf{2}(:,:,3)])));
mesh(1:nobs, freq_hat, real(spect_surf{2}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([minv, maxv]); title('$f_{22}$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(spect_surf{2}(:,:,2))); 
mesh(1:nobs, freq_hat, real(spect_surf{2}(:,:,3)));
subplot(1,3,3);
maxv = real(max(max([spect_surf{3}(:,:,1), spect_surf{3}(:,:,2), spect_surf{3}(:,:,3)])));
minv = real(min(min([spect_surf{3}(:,:,1), spect_surf{3}(:,:,2), spect_surf{3}(:,:,3)])));
mesh(1:nobs, freq_hat, real(spect_surf{3}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([minv, maxv]); title('$f_{33}$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(spect_surf{3}(:,:,2))); 
mesh(1:nobs, freq_hat, real(spect_surf{3}(:,:,3)));

figure
subplot(1,3,1);
mesh(1:nobs, freq_hat, real(coh_surf{1}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([0, 1]); title('$ \rho_{21}^2$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(coh_surf{1}(:,:,2))); 
mesh(1:nobs, freq_hat, real(coh_surf{1}(:,:,3)));
hold off
subplot(1,3,2);
mesh(1:nobs, freq_hat, real(coh_surf{2}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([0, 1]); title('$ \rho_{31}^2$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(coh_surf{2}(:,:,2))); 
mesh(1:nobs, freq_hat, real(coh_surf{2}(:,:,3)));
hold off
subplot(1,3,3);
mesh(1:nobs, freq_hat, real(coh_surf{3}(:,:,1)));set(gca,'YDir','reverse'); xlabel('time'); ylabel('freq'); zlim([0, 1]); title('$ \rho_{32}^2$','Interpreter','LaTex')
hold on
mesh(1:nobs, freq_hat, real(coh_surf{3}(:,:,2))); 
mesh(1:nobs, freq_hat, real(coh_surf{3}(:,:,3)));
hold off



