function [out, fitparams] = MultiSpect_diag(ts,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for the MCMC nonstationary multivariate spectrum analysis paper
% 07/02/2016
% Does the MCMC iterations for bivariate time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract information from the option parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If params is empty, then use the default paraemters
if  nargin==1
    params = OptsMultiSpect();
else
    params = varargin{1};
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','MATLAB:nearlySingularMatrix')
nloop = params.nloop;   %number of total MCMC iterations
nwarmup = params.nwarmup; %number of warmup period
nexp_max = params.nexp_max;   %Input maximum number of segments

global dimen nobs nbasis nBeta sigmasqalpha tau_up_limit ...
            prob_mm1 tmin var_inflate_1 var_inflate_2 options        
options = optimset('Display','off','GradObj','on','Hessian','on','MaxIter',10000,...
    'MaxFunEvals',10000,'TolFun',1e-5,'TolX',1e-5);

nbasis = params.nbasis;             %number of linear smoothing spline basis functions 
nBeta = nbasis + 1;    %number of coefficients for real part of components
sigmasqalpha = params.sigmasqalpha;       %smoothing parameter for alpha
tau_up_limit = params.tau_up_limit;       %the prior for smoothing parameters    
prob_mm1 = params.prob_mml;         %the probability of small jump, 1-prob of big jump
tmin = params.tmin;             %minimum number of observation in each segment
var_inflate_1 = params.v1(1);      
var_inflate_2 = params.v1(1);
nfreq_hat = params.nfreq;
freq_hat=(0:nfreq_hat)'/(2*nfreq_hat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = size(ts);       
nobs = dim(1);
dimen = dim(2);


spect_hat = cell(nexp_max,1);
for j=1:nexp_max
    spect_hat{j} = zeros(2,2,nfreq_hat+1,j,nloop+1);
end

nexp_curr = params.init;  %initialize the number of segments

%initialize tausq
for j=1:nexp_curr
	tausq_curr(:,j,1)=rand(4,1)*tau_up_limit;
end

%initilize the location of the changepoints for j=1:nexp_curr(1)
xi_curr = zeros(nexp_curr,1);
nseg_curr = zeros(nexp_curr,1);
for j=1:nexp_curr
    if nexp_curr==1
        xi_curr = nobs;
        nseg_curr = nobs;
    else
        if j==1
			nposs = nobs-nexp_curr*tmin+1;
			xi_curr(j) = tmin + unidrnd(nposs)-1;
			nseg_curr(j) = xi_curr(j);
        elseif j>1 && j<nexp_curr
			nposs=nobs-xi_curr(j-1)-tmin*(nexp_curr-j+1)+1;
			xi_curr(j)=tmin+unidrnd(nposs)+xi_curr(j-1)-1;
			nseg_curr(j)=xi_curr(j)-xi_curr(j-1);
        else
			xi_curr(j)=nobs;
			nseg_curr(j)=xi_curr(j)-xi_curr(j-1);	
        end
    end
end

%index matrix for which compoents of choloskey decomposition changed
chol_index = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1;...
              1 1 0 0; 1 0 1 0; 1 0 0 1; 0 0 1 1; 0 1 0 1;...
              0 1 1 0; 1 1 1 0; 1 1 0 1; 0 1 1 1; 1 0 1 1;...
              1 1 1 1];
% 1: no change; 2-5: one changes; 6-11: two changes; 12-15: three changes
% 16: all changed

Beta_curr = zeros(nBeta,4,nexp_curr);
Phi_temp = repmat(16,nexp_curr,1);
for j=1:nexp_curr
   	[Beta_mean, Beta_var,~] = postBeta1(chol_index, Phi_temp(j), j, ts, tausq_curr, Beta_curr, xi_curr);
	Beta_curr(:,:,j) = reshape(mvnrnd(Beta_mean,0.5*(Beta_var+Beta_var')),nBeta,4);
end
Phi_curr=[repmat(16,nexp_curr-1,1);0];
%jumping probabilities
epsilon=zeros(nloop,1);
met_rat=zeros(nloop,1);
rat_birth=[];
rat_death=[];


%%%initialize .mat file
S.('ts') = ts;
S.('varargin') = varargin;
S.('nexp') = [];
S.('nseg') = [];
S.('xi') = {};
S.('Beta') = {};
S.('tausq') = {};
S.('spect_hat') = {};
S.('Spec_1') = {};
S.('Spec_2') = {};
S.('Coh') = {};
S.('Phi') = {};
if params.convdiag==1
   S.('conv_diag') = []; 
end

save(params.fname, '-struct', 'S', '-v7.3');

%create matfile object to use in MCMC sampler
out = matfile(params.fname, 'Writable', true);

%create data structures to store temporary results in between ouputs
Beta_tmp = cell(params.batchsize,1);
spect_hat_tmp = cell(params.batchsize,1);
nexp_tmp = zeros(params.batchsize,1);
nseg_tmp = cell(params.batchsize,1);
tausq_tmp = cell(params.batchsize,1);
xi_tmp = cell(params.batchsize,1);
Phi_tmp = cell(params.batchsize,1);
Spect_1_tmp = cell(params.batchsize,1);
Spect_2_tmp = cell(params.batchsize,1);
Coh_tmp = cell(params.batchsize,1);
if params.convdiag==1
   convdiag_tmp = zeros(params.batchsize,4*nBeta+4); 
end    

% preallocate for the worst case memory use
for i=1:params.batchsize
    Beta_tmp{i} = zeros(nBeta, 4, params.nexp_max);
    spect_hat_tmp{i} = zeros(2,2,nfreq_hat+1,params.nexp_max);
    nexp_tmp(i) = params.nexp_max;
    nseg_tmp{i} = zeros(params.nexp_max,1);
    tausq_tmp{i} = zeros(4,params.nexp_max);
    xi_tmp{i} = zeros(params.nexp_max,1);
    Phi_tmp{i} = zeros(params.nexp_max,1);
    Spect_1_tmp{i} = zeros(nfreq_hat+1,nobs);
    Spect_2_tmp{i} = zeros(nfreq_hat+1,nobs);
    Coh_tmp{i} = zeros(nfreq_hat+1,nobs);
    if params.convdiag==1
        convdiag_tmp(i,:)=zeros(1,4*nBeta+4);
    end
end

%set batch index
batch_idx = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:nloop
    tic;
    if(mod(p,100)==0)
       fprintf('iter: %g of %g \n' ,p, nloop)
    end
    if p<nwarmup
        var_inflate_1=1;
        var_inflate_2=1;
    else
        var_inflate_1=1;
        var_inflate_2=1;
    end
    
    %========================
    %BETWEEN MODEL MOVE
    %========================
    
    kk = length(find(nseg_curr>2*tmin)); %Number of available segments
    
    %===========================
    %Deciding on birth or death
    if kk==0 %Stay where you (if nexp_curr=1) or join segments if there are no available segments to cut
        if nexp_curr==1
            nexp_prop = nexp_curr; %Stay 
            log_move_prop = 0;
            log_move_curr = 0;
        else
            nexp_prop = nexp_curr-1; %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 1;
            else
                log_move_curr = log(0.5);
            end
        end
    else
        if nexp_curr==1
            nexp_prop = nexp_curr + 1; %birth
            log_move_prop = 0;
            if nexp_prop==nexp_max
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        elseif nexp_curr==nexp_max
            nexp_prop = nexp_curr-1;   %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        else
            u = rand;
            if u<0.5
                nexp_prop = nexp_curr+1; %birth
                if nexp_prop==nexp_max;
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr=log(0.5);
                    log_move_prop=log(0.5);
                end
            else
                nexp_prop = nexp_curr-1; %death
                if nexp_prop==1
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr = log(0.5);
                    log_move_prop = log(0.5);
                end
            end
        end
    end

    
    if nexp_prop<nexp_curr
        %Do Death step
        [met_rat(p),nseg_prop,xi_prop,tausq_prop,Beta_prop, Phi_prop]= death(chol_index,ts,nexp_curr,nexp_prop,...
		tausq_curr,xi_curr,nseg_curr,Beta_curr,Phi_curr, log_move_curr,log_move_prop);
        rat_death=[rat_death met_rat(p)];
    elseif nexp_prop>nexp_curr
        %Do Birth step
        [met_rat(p),nseg_prop,xi_prop,tausq_prop,Beta_prop, Phi_prop]= birth(chol_index,ts,nexp_curr,nexp_prop,...
		tausq_curr,xi_curr,nseg_curr,Beta_curr,Phi_curr,log_move_curr,log_move_prop);
        rat_birth=[rat_birth met_rat(p)];
    else
        xi_prop=xi_curr;
        nseg_prop=nseg_curr;
        tausq_prop=tausq_curr;
        Beta_prop=Beta_curr;
        Phi_prop=Phi_curr;
        met_rat(p) = 1;
    end
    u = rand;
    if u<met_rat(p)
		nexp_curr=nexp_prop;
		xi_curr=xi_prop;
		nseg_curr=nseg_prop;
		tausq_curr=tausq_prop;
		Beta_curr=Beta_prop;
        Phi_curr=Phi_prop;
    end
    %========================
    %WITHIN MODEL MOVE
    %========================
	%Drawing a new cut point and Beta simultaneously
    %update coeffiecient of linear basis function

    [epsilon(p),nseg_new,xi_prop,~,Beta_prop,Phi_prop,seg_prop]= ...
    within(chol_index,ts,nexp_curr,tausq_curr, xi_curr,...
                                                            nseg_curr, Beta_curr, Phi_curr);
    u = rand;
    if (u<epsilon(p)|| p==1) 
        if nexp_curr>1
            for j=seg_prop:seg_prop+1
                Beta_curr=Beta_prop;
                xi_curr=xi_prop;
                nseg_curr=nseg_new;
                Phi_curr=Phi_prop;
            end
        else
            Beta_curr=Beta_prop;
        end
    end
    
    %Drawing tau
    for j=1:nexp_curr
        for i=1:3
            tau_a = nbasis/2;
            tau_b = sum(Beta_curr(2:end,i,j).^2)/2;
            u=rand;
            const1 = gamcdf(1/tau_up_limit,tau_a,1/tau_b);
            const2 = 1-u*(1-const1);
            tausq_curr(i,j) = 1/gaminv(const2,tau_a,1/tau_b);
        end
        tau_a = nbasis/2;
        tau_b = sum(Beta_curr(1:nBeta,4,j).^2)/2;
        u=rand;
        const1 = gamcdf(1/tau_up_limit,tau_a,1/tau_b);
        const2 = 1-u*(1-const1);
        tausq_curr(4,j) = 1/gaminv(const2,tau_a,1/tau_b);
    end
    
    %==================================
    %Estimating Spectral Density
    %==================================
    [xx_r, xx_i] = lin_basis_func(freq_hat); %produce linear basis functions
    
    for j =1:nexp_curr
        %getting the coefficients of linear basis functions
        g1 = Beta_curr(:,1:3,j);
        g2 = Beta_curr(1:nBeta,4,j);
        
        theta_real = xx_r * g1(:,3);
        theta_imag = xx_i * g2;
        theta = theta_real+sqrt(-1)*theta_imag;
        
        delta_sq_hat = zeros(2,nfreq_hat+1);
        for q=1:2
            delta_sq_hat(q,:) = exp(xx_r * g1(:,q))';
        end
        %produce the spectral density matrix
        for k=1:(nfreq_hat+1)
            TT = eye(2);
            TT(2,1) = -theta(k);
            spect_hat{nexp_curr}(:,:,k,j,p+1) = ...
            inv(TT)*diag(delta_sq_hat(:,k))*inv(TT');
        end
    end
    
    %============================================
    % Get time-varying spectra and coherence
    %============================================
    Spec_1_est = zeros(nfreq_hat+1,nobs);
    Spec_2_est = zeros(nfreq_hat+1,nobs);
    Coh_est = zeros(nfreq_hat+1,nobs);
    
    % first spectrum 
    if((p+1)>nwarmup)
        spec_hat_curr=squeeze(spect_hat{nexp_curr}(1,1,:,:,p+1));
        for j=1:nexp_curr
            if(j==1)
               Spect_1_est(:,1:xi_curr(j))= repmat(spec_hat_curr(:,j),1,xi_curr(j));
            else
               Spect_1_est(:,xi_curr(j-1)+1:xi_curr(j))= ...
                repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1));
            end
        end
    end
    % second spectrum 
    if((p+1)>nwarmup)
        spec_hat_curr=squeeze(spect_hat{nexp_curr}(2,2,:,:,p+1));
        for j=1:nexp_curr
            if(j==1)
               Spect_2_est(:,1:xi_curr(j))= repmat(spec_hat_curr(:,j),1,xi_curr(j));
            else
               Spect_2_est(:,xi_curr(j-1)+1:xi_curr(j))= ...
                repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1));
            end
        end
    end
    % coherence
    if((p+1)>nwarmup)
        spec_hat_curr=abs(squeeze(spect_hat{nexp_curr}(2,1,:,:,p+1))).^2./...
            (squeeze(spect_hat{nexp_curr}(1,1,:,:,p+1)).*squeeze(spect_hat{nexp_curr}(2,2,:,:,p+1)));
        for j=1:nexp_curr
            if(j==1)
               Coh_est(:,1:xi_curr(j))=repmat(spec_hat_curr(:,j),1,xi_curr(j));
            else
               Coh_est(:,xi_curr(j-1)+1:xi_curr(j))=repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1));
            end
        end
    end
    
    %create convergence diagnostics
    if params.convdiag==1
        convdiag_curr=convdiag(nobs, nBeta, nseg_curr, nexp_curr, Beta_curr, tausq_curr);
    end 
    
    %output to temporary container
    Beta_tmp{batch_idx} = Beta_curr;
    spect_hat_tmp{batch_idx} = spect_hat;
    nexp_tmp(batch_idx) = nexp_curr;
    nseg_tmp{batch_idx} = nseg_curr;
    tausq_tmp{batch_idx} = tausq_curr;
    xi_tmp{batch_idx} = xi_curr;
    Phi_tmp{batch_idx} = Phi_curr;
    if((p+1)>nwarmup)
        Spect_1_tmp{batch_idx} = Spect_1_est;
        Spect_2_tmp{batch_idx} = Spect_2_est;
        Coh_tmp{batch_idx} = Coh_est;
    end    
    if params.convdiag==1
        convdiag_tmp=convdiag_curr;
    end   
    
    if(mod(p,params.batchsize)==0)
        
        
        %output data to .mat file
        out.Beta((p-params.batchsize+1):p,1) = Beta_tmp;
        out.spect_hat((p-params.batchsize+1):p,1) = spect_hat_tmp;    
        out.nexp((p-params.batchsize+1):p,1)=nexp_tmp;
        out.nseg((p-params.batchsize+1):p,1)=nseg_tmp;
        out.tausq((p-params.batchsize+1):p,1) = tausq_tmp;
        out.xi((p-params.batchsize+1):p,1) = xi_tmp;
        out.Phi((p-params.batchsize+1):p,1) = Phi_tmp;
        out.Spect_1((p-params.batchsize+1):p,1) = Spect_1_tmp;    
        out.Spect_2((p-params.batchsize+1):p,1) = Spect_2_tmp; 
        out.Coh((p-params.batchsize+1):p,1) = Coh_tmp;        
        if params.convdiag==1  
            out.convdiag((p-params.batchsize+1):p,1:nbeta+1)=convdiag_tmp;
        end  
        
        %reset temporary data containers        
        Beta_tmp = cell(params.batchsize,1);
        spect_hat_tmp = cell(params.batchsize,1);
        nexp_tmp = zeros(params.batchsize,1);
        nseg_tmp = cell(params.batchsize,1);
        tausq_tmp = cell(params.batchsize,1);
        xi_tmp = cell(params.batchsize,1);
        Phi_tmp = cell(params.batchsize,1);
        Spect_1_tmp = cell(params.batchsize,1);
        Spect_2_tmp = cell(params.batchsize,1);
        Coh_tmp = cell(params.batchsize,1);
        if params.convdiag==1
            convdiag_tmp = zeros(params.batchsize,4*nBeta+4); 
        end 
       
        %reset batch index
        batch_idx=0;
        
    end
    
    %increment batch index
    batch_idx = batch_idx + 1;
    
    tms(p) = toc;
    if params.verb ==1
        fprintf('sec / min hr: %g %g %g \n' ,[tms(p),sum(tms(1:p))/60, sum(tms(1:p))/(60*60)]')          
    end
    
end
fitparams = struct('nloop', nloop, ...
                   'nwarmup', nwarmup, ...
                   'timeMean', mean(tms), ...
                   'timeMax', max(tms(2:end)), ...
                   'timeMin', min(tms),...
                   'timeStd', std(tms(2:end)));
end



