%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for the MCMC nonstationary multivariate spectrum analysis paper
% 07/02/2016
% Does the MCMC iterations for bivariate time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   I) get the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  load('abrupt_example.mat') 
load('abrupt_example2.mat') 

figure
subplot(2,1,1);
plot(zt(:,1),'color','black')
subplot(2,1,2);
plot(zt(:,1),'color','black')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','MATLAB:nearlySingularMatrix')
nloop = 6000;   %number of total MCMC iterations
nwarmup = 2000; %number of warmup period
nexp_max = 5;   %Input maximum number of segments

global dimen nobs nbasis ngamma sigmasqalpha tau_up_limit ...
            prob_mm1 tmin var_inflate_1 var_inflate_2 options        
options = optimset('Display','off','GradObj','on','Hessian','on','MaxIter',10000,...
    'MaxFunEvals',10000,'TolFun',1e-5,'TolX',1e-5);

ts = zt;               %a T by two time series 
dim = size(ts);       
nobs = dim(1);
dimen = dim(2);
n = nobs;
yy = fft(ts)/sqrt(n);  %DFT 
y = yy(1:ceil(n/2),:);
nf = length(y);
tt=linspace(1,nobs,nobs)';

nbasis = 7;             %number of linear smoothing spline basis functions 
ngamma = nbasis + 1;    %number of coefficients for real part of components
sigmasqalpha = 5;       %smoothing parameter for alpha
tau_up_limit = 5;       %the prior for smoothing parameters    
prob_mm1 = 0.8;         %the probability of small jump, 1-prob of big jump
tmin = 100;             %minimum number of observation in each segment
var_inflate_1 = 2.5;      
var_inflate_2 = 2.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% II) Run the estimation proceedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfreq_hat = 50;
freq_hat=(0:nfreq_hat)'/(2*nfreq_hat);

tausq = cell(nexp_max,1);  %big array for smoothing parameters
gamma = cell(nexp_max,1);  %big array for coefficients 
spect_hat = cell(nexp_max,1);
xi = cell(nexp_max,1);     %Cutpoint locations xi_1 is first cutpoint, xi_) is beginning of timeseries
nseg = cell(nexp_max,1);   %Number of observations in each segment
rho = cell(nexp_max,1);    %index for component change
spec_hat = cell(nexp_max,1);
fmean = cell(nexp_max,1);

for j=1:nexp_max
    tausq{j}=ones(4,j,nloop+1);
    gamma{j} = zeros(ngamma,4,j,nloop+1);
    spect_hat{j} = zeros(2,2,nfreq_hat+1,j,nloop+1);
    fmean{j} = zeros(dimen,dimen,nfreq_hat+1,j);
    xi{j}=ones(j,nloop+1);
	nseg{j}=ones(j,nloop+1);
    rho{j}=zeros(j,nloop+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilize the MCMC iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexp_curr=nexp_max*ones(nloop+1,1); %big array for number of segment
nexp_curr(1) = 3;  %initialize the number of segments

%initialize tausq
for j=1:nexp_curr(1)
	tausq{nexp_curr(1)}(:,j,1)=rand(4,1)*tau_up_limit;
end

%initilize the location of the changepoints for j=1:nexp_curr(1)
for j=1:nexp_curr(1)
    if nexp_curr(1)==1
        xi{nexp_curr(1)}(j,1) = nobs;
        nseg{nexp_curr(1)}(j,1) = nobs;
    else
        if j==1
			nposs = nobs-nexp_curr(1)*tmin+1;
			xi{nexp_curr(1)}(j,1) = tmin + unidrnd(nposs)-1;
			nseg{nexp_curr(1)}(j,1) = xi{nexp_curr(1)}(j,1);
        elseif j>1 && j<nexp_curr(1)
			nposs=nobs-xi{nexp_curr(1)}(j-1,1)-tmin*(nexp_curr(1)-j+1)+1;
			xi{nexp_curr(1)}(j,1)=tmin+unidrnd(nposs)+xi{nexp_curr(1)}(j-1,1)-1;
			nseg{nexp_curr(1)}(j,1)=xi{nexp_curr(1)}(j,1)-xi{nexp_curr(1)}(j-1,1);
        else
			xi{nexp_curr(1)}(j,1)=nobs;
			nseg{nexp_curr(1)}(j,1)=xi{nexp_curr(1)}(j,1)-xi{nexp_curr(1)}(j-1,1);	
        end
    end
end

%initilize indicator for order of rho occur.
Z = zeros(nobs,1);
Z(xi{nexp_curr(1)}(:,1)) = Z(xi{nexp_curr(1)}(:,1)) + 1;

%index matrix for which compoents of choloskey decomposition changed
chol_index = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1;...
              1 1 0 0; 1 0 1 0; 1 0 0 1; 0 0 1 1; 0 1 0 1;...
              0 1 1 0; 1 1 1 0; 1 1 0 1; 0 1 1 1; 1 0 1 1;...
              1 1 1 1];
% 1: no change; 2-5: one changes; 6-11: two changes; 12-15: three changes
% 16: all changed

xi_temp = xi{nexp_curr(1)}(:,1);
nseg_temp = nseg{nexp_curr(1)}(:,1);
tau_temp = tausq{nexp_curr(1)}(:,:,1);
gamma_temp = gamma{nexp_curr(1)}(:,:,:,1);
rho{nexp_curr(1)}(1:(end-1),1)=16;
rho_temp = repmat(16,nexp_curr(1));
for j=1:nexp_curr(1)
   	[gamma_mean, gamma_var,~] = postgamma1(chol_index, rho_temp(j), j, ts, tau_temp(:,j), gamma_temp(:,:,j,1), xi_temp);
	gamma{nexp_curr(1)}(:,:,j,1) = reshape([mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')),0],ngamma,4);
end
%jumping probabilities
epsilon=zeros(nloop,1);
met_rat=zeros(nloop,1);

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:nloop
    if(mod(p,100)==0)
        p
    end
    
    if p<nwarmup
        var_inflate_1=2.5;
        var_inflate_2=2.5;
    else
        var_inflate_1=2;
        var_inflate_2=2;
    end
    
    %========================
    %BETWEEN MODEL MOVE
    %========================
    
    kk = length(find(nseg{nexp_curr(p)}(:,p)>2*tmin)); %Number of available segments
    
    %===========================
    %Deciding on birth or death
    if kk==0 %Stay where you (if nexp_curr=1) or join segments if there are no available segments to cut
        if nexp_curr(p)==1
            nexp_prop = nexp_curr(p); %Stay 
            log_move_prop = 0;
            log_move_curr = 0;
        else
            nexp_prop = nexp_curr(p)-1; %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 1;
            else
                log_move_curr = log(0.5);
            end
        end
    else
        if nexp_curr(p)==1
            nexp_prop = nexp_curr(p) + 1; %birth
            log_move_prop = 0;
            if nexp_prop==nexp_max
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        elseif nexp_curr(p)==nexp_max
            nexp_prop = nexp_curr(p)-1;   %death
            log_move_prop = 0;
            if nexp_prop==1
                log_move_curr = 0;
            else
                log_move_curr = log(0.5);
            end
        else
            u = rand;
            if u<0.5
                nexp_prop = nexp_curr(p)+1; %birth
                if nexp_prop==nexp_max;
                    log_move_curr = 0;
                    log_move_prop = log(0.5);
                else
                    log_move_curr=log(0.5);
                    log_move_prop=log(0.5);
                end
            else
                nexp_prop = nexp_curr(p)-1; %death
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
    xi_curr_temp = xi{nexp_curr(p)}(:,p);
    gamma_curr_temp = gamma{nexp_curr(p)}(:,:,:,p);
    nseg_curr_temp = nseg{nexp_curr(p)}(:,p);
    tau_curr_temp = tausq{nexp_curr(p)}(:,:,p);
    rho_curr_temp = rho{nexp_curr(p)}(:,p);
    Z(xi{nexp_curr(p)}(:,p)) = Z(xi{nexp_curr(p)}(:,p)) + 1;
    
    if nexp_prop<nexp_curr(p)
        %Do Death step
        [met_rat(p),nseg_prop,xi_prop,tau_prop,gamma_prop, rho_prop, z_tmp]= death(chol_index,ts,nexp_curr(p),nexp_prop,...
		tau_curr_temp,xi_curr_temp,nseg_curr_temp,gamma_curr_temp,rho_curr_temp, log_move_curr,log_move_prop,Z);
    elseif nexp_prop>nexp_curr(p)
        %Do Birth step
        [met_rat(p),nseg_prop,xi_prop,tau_prop,gamma_prop, rho_prop, z_tmp]= birth(chol_index,ts,nexp_curr(p),nexp_prop,...
		tau_curr_temp,xi_curr_temp,nseg_curr_temp,gamma_curr_temp,rho_curr_temp,log_move_curr,log_move_prop,Z);
    else
        xi_prop=xi{nexp_curr(p)}(:,p);
        nseg_prop=nseg{nexp_curr(p)}(:,p);
        tau_prop=tausq{nexp_curr(p)}(:,:,p);
        gamma_prop=beta{nexp_curr(p)}(:,:,p);
        rho_prop=rho{nexp_curr(p)}(:,p);
        met_rat(p) = 1;
    end
    u = rand;
    if u<met_rat(p)
		nexp_curr(p+1)=nexp_prop;
		xi{nexp_curr(p+1)}(:,p+1)=xi_prop;
		nseg{nexp_curr(p+1)}(:,p+1)=nseg_prop;
		tausq{nexp_curr(p+1)}(:,:,p+1)=tau_prop;
		gamma{nexp_curr(p+1)}(:,:,:,p+1)=gamma_prop;
        rho{nexp_curr(p+1)}(:,p+1)=rho_prop;
        Z=z_tmp;
    else
		nexp_curr(p+1)=nexp_curr(p);
		xi{nexp_curr(p+1)}(:,p+1)=xi{nexp_curr(p+1)}(:,p);
		nseg{nexp_curr(p+1)}(:,p+1)=nseg{nexp_curr(p+1)}(:,p);
		tausq{nexp_curr(p+1)}(:,:,p+1)=tausq{nexp_curr(p+1)}(:,:,p);
		gamma{nexp_curr(p+1)}(:,:,:,p+1)=gamma{nexp_curr(p+1)}(:,:,:,p);
        rho{nexp_curr(p+1)}(:,p+1)=rho{nexp_curr(p+1)}(:,p);
    end
    %========================
    %WITHIN MODEL MOVE
    %========================
	%Drawing a new cut point and gamma simultaneously
    %update coeffiecient of linear basis function
    xi_curr_temp=xi{nexp_curr(p+1)}(:,p+1);
    gamma_curr_temp=gamma{nexp_curr(p+1)}(:,:,:,p+1);
    tau_temp=tausq{nexp_curr(p+1)}(:,:,p+1);
    nseg_curr_temp=nseg{nexp_curr(p+1)}(:,p+1);
    rho_temp = rho{nexp_curr(p+1)}(:,p+1); 
    [epsilon(p),nseg_new,xi_prop,tau_prop,gamma_prop,rho_prop,seg_temp,z_tmp]= ...
    within(chol_index,ts,nexp_curr(p+1),tau_temp, xi_curr_temp,...
                                                            nseg_curr_temp, gamma_curr_temp,rho_temp,Z);
    u = rand;
    if (u<epsilon(p)|| p==1) 
        if nexp_curr(p+1)>1
            for j=seg_temp:seg_temp+1
                gamma{nexp_curr(p+1)}(:,:,j,p+1)=gamma_prop(:,:,j);
                xi{nexp_curr(p+1)}(j,p+1)=xi_prop(j);
                nseg{nexp_curr(p+1)}(j,p+1)=nseg_new(j);
                rho{nexp_curr(p+1)}(j,p+1)=rho_prop(j);
                Z=z_tmp;
            end
        else
            gamma{nexp_curr(p+1)}(:,1,p+1)=gamma_prop(:,1);
        end
    else
        gamma{nexp_curr(p+1)}(:,:,:,p+1)=gamma_curr_temp;
        xi{nexp_curr(p+1)}(:,p+1)=xi_curr_temp;
        nseg{nexp_curr(p+1)}(:,p+1)=nseg_curr_temp;
        rho{nexp_curr(p+1)}(:,p+1)=rho_temp;
    end
    
    %Drawing tau
    for j=1:nexp_curr(p+1)
        for i=1:3
            tau_a = nbasis/2;
            tau_b = sum(gamma{nexp_curr(p+1)}(2:end,i,j,p+1).^2)/2;
            u=rand;
            const1 = gamcdf(1/tau_up_limit,tau_a,1/tau_b);
            const2 = 1-u*(1-const1);
            tausq{nexp_curr(p+1)}(i,j,p+1) = 1/gaminv(const2,tau_a,1/tau_b);
%           tausq{nexp_curr(p+1)}(i,j,p+1)=1/gamrnd(tau_a,1/tau_b);
        end
        tau_a = nbasis/2;
        tau_b = sum(gamma{nexp_curr(p+1)}(1:nbasis,4,j,p+1).^2)/2;
        u=rand;
        const1 = gamcdf(1/tau_up_limit,tau_a,1/tau_b);
        const2 = 1-u*(1-const1);
        tausq{nexp_curr(p+1)}(4,j,p+1) = 1/gaminv(const2,tau_a,1/tau_b);
%       tausq{nexp_curr(p+1)}(4,j,p+1)=1/gamrnd(tau_a,1/tau_b);
    end
    
    %==================================
    %Estimating Spectral Density
    %==================================
    [xx_r, xx_i] = lin_basis_func(freq_hat); %produce linear basis functions
    
    for j =1:nexp_curr(p+1)
        %getting the coefficients of linear basis functions
        g1 = gamma{nexp_curr(p+1)}(:,1:3,j,p+1);
        g2 = gamma{nexp_curr(p+1)}(1:nbasis,4,j,p+1);
        
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
            spect_hat{nexp_curr(p+1)}(:,:,k,j,p+1) = ...
            inv(TT)*diag(delta_sq_hat(:,k))*inv(TT');
        end
    end    
end

run_time = toc/60;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average ove MCMC samples to get time-varying spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot the estimated of time-varying spectrum (f11)
s_11=zeros(nfreq_hat+1,nobs);
for p=1:nloop
    if(p>nwarmup)
        xi_curr=xi{nexp_curr(p)}(:,p);
        spec_hat_curr=squeeze(spect_hat{nexp_curr(p)}(1,1,:,:,p));
        for j=1:nexp_curr(p)
            if(j==1)
               s_11(:,1:xi_curr(j))=s_11(:,1:xi_curr(j))+repmat(spec_hat_curr(:,j),1,xi_curr(j))/(nloop-nwarmup);
            else
               s_11(:,xi_curr(j-1)+1:xi_curr(j))=s_11(:,xi_curr(j-1)+1:xi_curr(j))+...
                repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
            end
        end
    end
end
figure
h_11=surf(tt,freq_hat,real(s_11));
set(h_11, 'edgecolor','none')
set(gca,'YDir','reverse')

%plot the estimated of time-varying spectrum (f22)
s_22=zeros(nfreq_hat+1,nobs);
for p=1:nloop
    if(p>nwarmup)
        xi_curr=xi{nexp_curr(p)}(:,p);
        spec_hat_curr=squeeze(spect_hat{nexp_curr(p)}(2,2,:,:,p));
        for j=1:nexp_curr(p)
            if(j==1)
               s_22(:,1:xi_curr(j))=s_22(:,1:xi_curr(j))+repmat(spec_hat_curr(:,j),1,xi_curr(j))/(nloop-nwarmup);
            else
               s_22(:,xi_curr(j-1)+1:xi_curr(j))=s_22(:,xi_curr(j-1)+1:xi_curr(j))+...
                repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
            end
        end
    end
end
figure
h_22=surf(tt,freq_hat,real(s_22));
set(h_22, 'edgecolor','none')
set(gca,'YDir','reverse')

%plot the estimated of time-varying spectrum (f21)
s_21=zeros(nfreq_hat+1,nobs);
for p=1:nloop
    if(p>nwarmup)
        xi_curr=xi{nexp_curr(p)}(:,p);
        spec_hat_curr=squeeze(spect_hat{nexp_curr(p)}(2,1,:,:,p));
        for j=1:nexp_curr(p)
            if(j==1)
               s_21(:,1:xi_curr(j))=s_21(:,1:xi_curr(j))+repmat(spec_hat_curr(:,j),1,xi_curr(j))/(nloop-nwarmup);
            else
               s_21(:,xi_curr(j-1)+1:xi_curr(j))=s_21(:,xi_curr(j-1)+1:xi_curr(j))+...
                repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
            end
        end
    end
end
s_21_real = real(s_21);
s_21_imag = imag(s_21);

figure
h_21_real=surf(tt,freq_hat,real(s_21_real));
set(h_21_real, 'edgecolor','none')
set(gca,'YDir','reverse')

figure
h_21_imag=surf(tt,freq_hat,real(s_21_imag));
set(h_21_imag, 'edgecolor','none')
set(gca,'YDir','reverse')

%plot the estimated of coherence
s_coh=zeros(nfreq_hat+1,nobs);
for p=1:nloop
    if(p>nwarmup)
        xi_curr=xi{nexp_curr(p)}(:,p);
        spec_hat_curr=abs(squeeze(spect_hat{nexp_curr(p)}(2,1,:,:,p))).^2./...
            (squeeze(spect_hat{nexp_curr(p)}(1,1,:,:,p)).*squeeze(spect_hat{nexp_curr(p)}(2,2,:,:,p)));
        for j=1:nexp_curr(p)
            if(j==1)
               s_coh(:,1:xi_curr(j))=s_coh(:,1:xi_curr(j))+repmat(spec_hat_curr(:,j),1,xi_curr(j))/(nloop-nwarmup);
            else
               s_coh(:,xi_curr(j-1)+1:xi_curr(j))=s_coh(:,xi_curr(j-1)+1:xi_curr(j))+...
                repmat(spec_hat_curr(:,j),1,xi_curr(j)-xi_curr(j-1))/(nloop-nwarmup);
            end
        end
    end
end
figure
h_coh=surf(tt,freq_hat,real(s_coh));
set(h_coh, 'edgecolor','none')
set(gca,'YDir','reverse')

%partitions
figure
histogram(nexp_curr(nwarmup+1:nloop),'FaceAlpha',1, 'Normalization','probability','BinMethod','integers','BinWidth',.5,...
    'BinLimits',[1,10])
for j=1:nexp_max
	kk=find(nexp_curr(nwarmup+1:nloop)==j);
    post_prob_nexp(j)=length(kk)/(nloop-nwarmup);
    if ~isempty(kk) && j>1
        figure
        hold
        title(['Plot of Partition Points for a mixture of ',int2str(j)])
        for k=1:j-1
            plot(xi{j}(k,kk+nwarmup))
        end
        for k=1:j-1
            figure
            hold
            title(['Histogram for a mixture of ',int2str(j)])
            histogram(xi{j}(k,kk+nwarmup),'FaceAlpha',1, 'Normalization','probability','BinMethod','integers','BinWidth',20,...
    'BinLimits',[1,nobs])
        end
    end
%     if j==2
%         c1=min(xi{j}(1,kk+nwarmup));
%         c2=max(xi{j}(1,kk+nwarmup));
%         cut_grid=linspace(c1,c2,c2-c1)';
%         sort_cut=sort(xi{j}(1,kk+nwarmup));
%     [nout xout]= hist(sort_cut,cut_grid);
%     pexp=cumsum(nout./sum(nout));
%     figure('Visible','Off')
%     plot(xout,pexp)
end


