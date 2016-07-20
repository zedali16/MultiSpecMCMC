function [f,gr,h] = gamma_derive1(x, yobs_tmp, chol_index, rho_temp, tau_temp, gamma_temp, sigmasqalpha, nbasis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/04/2016 for the MCMC Non-Multi-spectrum analysis
% Function used for optimization process for those coeeficient selected to
% be different
%
%   Input:
%       1) x - initial values for coeeficient of basis functions need to be
%       optimized 
%       2) yobs_tmp - time series data within the segment
%       3) chol_index - index matrix
%       4) rho_temp - which component changed
%       5) tau_temp - smoothing parameters  
%       6) gamma_temp - current coefficients
%       7) sigmasqalpha - smoothing parameters for the constant in real
%       components
%       8) nbasis - number of basis function used
%   Main Outputs:
%       1) f - log posterior probability based on input parameters
%       2) gr - gradients for optimization process
%       3) h - hession matrix for optimization process
%
%   Required programs: lin_basis_func

%initilize gamma_1 and gamma_2: gamma_2 is for imaginary components
ngamma = nbasis + 1;
gamma_1 = zeros(ngamma,3);
gamma_2 = zeros(nbasis,1);

select = chol_index(rho_temp,:).*[1 2 3 4];
select = select(select~=0);

if isempty(find(select==4,1))~=0
   x = reshape(x,ngamma,length(select));
   gamma_temp(:,select) = x;
else
   x1 = reshape(x(1:((length(select)-1)*ngamma)), ngamma, (length(select)-1));
   x2 = reshape(x(((length(select)-1)*ngamma+1):end), nbasis, 1);
   gamma_temp(:,select(1:(end-1))) = x1;
   gamma_temp(1:nbasis,4) = x2;
end

% gamma_temp(:,select) = 0;
gamma_1(1:ngamma,1:3) = gamma_temp(1:ngamma,1:3);
gamma_2(1:nbasis) = gamma_temp(1:nbasis,4);

dim = size(yobs_tmp);
n = dim(1);
nfreq = floor(n/2);
tt = (0:nfreq)/(2*nfreq);
yy = fft(yobs_tmp)/sqrt(n);
y = yy(1:(nfreq+1),:);
nf = length(y);

delta_sq = zeros(2,nf);
[xx_r, xx_i]=lin_basis_func(tt);

%theta's
theta_real = xx_r * gamma_1(:,3);
theta_imag = xx_i * gamma_2;
theta = theta_real + sqrt(-1)*theta_imag;

%delta's
delta_sq(1,:) = exp(xx_r * gamma_1(:,1));
delta_sq(2,:) = exp(xx_r * gamma_1(:,2));

if (mod(n,2)==1) %odd n
    f = -sum(log(delta_sq(1,2:end))' + log(delta_sq(2,2:end))' + abs(y(2:end,1)).^2.*exp(-xx_r(2:end,:)*gamma_1(:,1)) + ...
           abs(y(2:end,2) - theta(2:end).*y(2:end,1)).^2.*exp(-xx_r(2:end,:)*gamma_1(:,2))) - ...
           0.5*(log(delta_sq(1,1))' + log(delta_sq(2,1))' + abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1)) + ...
           abs(y(1,2) - theta(1).*y(1,1))^2.*exp(-xx_r(1,:)*gamma_1(:,2)));
else
    f = -sum(log(delta_sq(1,2:nfreq))' + log(delta_sq(2,2:nfreq))' + abs(y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*gamma_1(:,1)) + ...
           abs(y(2:nfreq,2) - theta(2:nfreq).*y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*gamma_1(:,2))) - ...
           0.5*(log(delta_sq(1,1)) + log(delta_sq(2,1)) + abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1)) + ...
           abs(y(1,2) - theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))) - ...
           0.5*(log(delta_sq(1,end)) + log(delta_sq(2,end)) + abs(y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,1)) + ...
           abs(y(end,2) - theta(end).*y(end,1))^2.*exp(-xx_r(end,:)*gamma_1(:,2)));
end
       
f = f - (0.5.*(gamma_1(1,1)* gamma_1(1,1)')/sigmasqalpha + 0.5.*(gamma_1(2:ngamma,1)'*gamma_1(2:ngamma,1))/tau_temp(1))*chol_index(rho_temp,1) -...
        (0.5.*(gamma_1(1,2)* gamma_1(1,2)')/sigmasqalpha + 0.5.*(gamma_1(2:ngamma,2)'*gamma_1(2:ngamma,2))/tau_temp(2))*chol_index(rho_temp,2) -...
        (0.5.*(gamma_1(1,3)* gamma_1(1,3)')/sigmasqalpha + 0.5.*(gamma_1(2:ngamma,3)'*gamma_1(2:ngamma,3))/tau_temp(3))*chol_index(rho_temp,3) -...
        (0.5.*(gamma_2(1:nbasis,1)'*gamma_2(1:nbasis,1))/tau_temp(4))*chol_index(rho_temp,4);

gr1 = zeros(ngamma,1);
gr2 = zeros(ngamma,1);
gr3 = zeros(ngamma,1);
gr4 = zeros(nbasis,1);

gr1(1) = gamma_1(1,1)/sigmasqalpha;
gr1(2:ngamma,1) = gamma_1(2:ngamma,1)/tau_temp(1);
gr2(1) = gamma_1(1,2)/sigmasqalpha;
gr2(2:ngamma,1) = gamma_1(2:ngamma,2)/tau_temp(2);
gr3(1) = gamma_1(1,3)/sigmasqalpha;
gr3(2:ngamma,1) = gamma_1(2:ngamma,3)/tau_temp(3);
gr4(1:nbasis,1) = gamma_2(1:nbasis,1)/tau_temp(4);

h11(1,1)=1/sigmasqalpha;
h11(2:ngamma,2:ngamma)=1/tau_temp(1)*eye(nbasis);
h22(1,1)=1/sigmasqalpha;
h22(2:ngamma,2:ngamma)= 1/tau_temp(2)*eye(nbasis);
h33(1,1)=1/sigmasqalpha;
h33(2:ngamma,2:ngamma)= 1/tau_temp(3)*eye(nbasis);
h44(1:nbasis,1:nbasis)=1/tau_temp(4)*eye(nbasis);
h42 = zeros(nbasis,ngamma);
h32 = zeros(ngamma,ngamma);

if (mod(n,2)==1)     
    %%%%%%%%%%%%%%%%%%%%%%%%
    %gradient
    %%%%%%%%%%%%%%%%%%%%%%%%
    rk = -y(2:end,1).*conj(y(2:end,2)) - y(2:end,2).*conj(y(2:end,1));
    ik = sqrt(-1)*(-y(2:end,1).*conj(y(2:end,2)) + y(2:end,2).*conj(y(2:end,1)));
    ck = 2*abs(y(2:end,1)).^2;
    
    gr1 = gr1 + xx_r(2:end,:)'*(1-abs(y(2:end,1)).^2.*exp(-xx_r(2:end,:)*gamma_1(:,1))) + ...
              0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))));
    gr2 = gr2 + xx_r(2:end,:)'*(1 - abs(y(2:end,2)-theta(2:end).*y(2:end,1)).^2.*exp(-xx_r(2:end,:)*gamma_1(:,2))) + ...
              0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))));  
    temp_mat_31 = bsxfun(@times, xx_r(2:end,:),rk);
    temp_mat_32 = bsxfun(@times, ck, bsxfun(@times, xx_r(2:end,:),xx_r(2:end,:)*gamma_1(:,3)));
    gr3 = gr3 + sum( bsxfun(@times, (temp_mat_31+temp_mat_32), exp(-xx_r(2:end,:)*gamma_1(:,2)) ))' +...
          0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:)' +...
          exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)');
    temp_mat_41 = bsxfun(@times, ik, xx_i(2:end,:));
    temp_mat_42 = bsxfun(@times, ck, bsxfun(@times,xx_i(2:end,:),xx_i(2:end,:)*gamma_2(:,1)));
    gr4 = gr4 + sum( bsxfun(@times, (temp_mat_41 + temp_mat_42), exp(-xx_r(2:end,:)*gamma_1(:,2))))' + ...
          0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1))))*xx_i(1,:)' +...
          exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)');
      
    %%%%%%%%%%%%%%%%%%%
    %hession
    %%%%%%%%%%%%%%%%%%%
    bigmat_h11 = kron(bsxfun(@times, abs(y(2:end,1)).^2.*exp(-xx_r(2:end,:)*gamma_1(:,1)), xx_r(2:end,:)),...
                                                                                ones(ngamma,1)');
    coefmat_h11 = repmat(xx_r(2:end,:), 1,ngamma);
    h11 = h11 + reshape(sum(bsxfun(@times, bigmat_h11, coefmat_h11),1),ngamma,ngamma) +...
            0.5*(abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))*xx_r(1,:)'*xx_r(1,:));
    
    bigmat_h22 = kron(bsxfun(@times,  abs(y(2:end,2)-theta(2:end).*y(2:end,1)).^2.*exp(-xx_r(2:end,:)*gamma_1(:,2)),...
                                                                xx_r(2:end,:)), ones(ngamma,1)');
    coefmat_h22 = repmat(xx_r(2:end,:), 1,ngamma);
    h22 = h22 + reshape(sum(bsxfun(@times, bigmat_h22, coefmat_h22),1),ngamma,ngamma) +...
            0.5*(abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))*xx_r(1,:)'*xx_r(1,:));
    
    bigmat_h33 = kron(bsxfun(@times, exp(-xx_r(2:end,:)*gamma_1(:,2)).*ck, xx_r(2:end,:)),...
                                                                                ones(ngamma,1)');
    coefmat_h33 = repmat(xx_r(2:end,:), 1,ngamma);         
    h33 = h33 + reshape(sum(bsxfun(@times, bigmat_h33, coefmat_h33),1),ngamma,ngamma) +...
                0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*xx_r(1,:)'*xx_r(1,:));


    bigmat_h44 = kron(bsxfun(@times, exp(-xx_r(2:end,:)*gamma_1(:,2)).*ck, xx_i(2:end,:)),...
                                                                                ones(nbasis,1)');
    coefmat_h44 = repmat(xx_i(2:end,:), 1,nbasis);         
    h44 = h44 + reshape(sum(bsxfun(@times, bigmat_h44, coefmat_h44),1),nbasis,nbasis) +...
                0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*xx_i(1,:)'*xx_i(1,:));
    
    bigmat_h42_1 = kron(bsxfun(@times, exp(-xx_r(2:end,:)*gamma_1(:,2)).*ck.*(xx_i(2:end,:)*gamma_2(:,1)),...
                        xx_i(2:end,:)), ones(ngamma,1)');
    coefmat_h42_1 = repmat(xx_r(2:end,:), 1,nbasis);  
    bigmat_h42_2 = kron(bsxfun(@times,exp(-xx_r(2:end,:)*gamma_1(:,2)).*ik, xx_i(2:end,:)), ones(ngamma,1)');
    coefmat_h42_2 = repmat(xx_r(2:end,:), 1,nbasis);  
    h42 = h42 + reshape(sum(bsxfun(@times, bigmat_h42_1, coefmat_h42_1) +...
                        bsxfun(@times, bigmat_h42_2, coefmat_h42_2),1),ngamma,nbasis)' +...
                        0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(2*abs(y(1,1)).^2*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)+...
                        sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1)))*xx_i(1,:))'*xx_r(1,:));

    bigmat_h32_1 = kron(bsxfun(@times, exp(-xx_r(2:end,:)*gamma_1(:,2)).*ck.*(xx_r(2:end,:)*gamma_1(:,3)),...
                        xx_r(2:end,:)), ones(ngamma,1)');
    coefmat_h32_1 = repmat(xx_r(2:end,:), 1,ngamma);  
    bigmat_h32_2 = kron(bsxfun(@times,exp(-xx_r(2:end,:)*gamma_1(:,2)).*rk, xx_r(2:end,:)), ones(ngamma,1)');
    coefmat_h32_2 = repmat(xx_r(2:end,:), 1,ngamma);  
    h32 = h32 + reshape(sum(bsxfun(@times, bigmat_h32_1, coefmat_h32_1) +...
                        bsxfun(@times, bigmat_h32_2, coefmat_h32_2),1),ngamma,ngamma)' +...
                        0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(2*abs(y(1,1)).^2*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)+...
                        (-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:))'*xx_r(1,:));      
%     for k=2:(nfreq+1)
%         
%         rk = -y(k,1).*conj(y(k,2)) - y(k,2).*conj(y(k,1));
%         ik = sqrt(-1)*(-y(k,1).*conj(y(k,2)) + y(k,2).*conj(y(k,1)));
%         ck = 2*abs(y(k,1)).^2;
%         gr1 = gr1 + xx_r(k,:)'*(1-abs(y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,1)));
%         gr2 = gr2 + xx_r(k,:)'*(1 - abs(y(k,2)-theta(k).*y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,2)));  
%         gr3 = gr3 + exp(-xx_r(k,:)*gamma_1(:,2))*rk*xx_r(k,:)' +...
%                     exp(-xx_r(k,:)*gamma_1(:,2))*ck*(xx_r(k,:)*gamma_1(:,3))*xx_r(k,:)';
%         gr4 = gr4 + exp(-xx_r(k,:)*gamma_1(:,2))*ik*xx_i(k,:)' +...
%                     exp(-xx_r(k,:)*gamma_1(:,2))*ck*(xx_i(k,:)*gamma_2(:,1))*xx_i(k,:)';
% 
%         %%%%%%%%%%%%%%%%%%%%
%         %hession        
%         %%%%%%%%%%%%%%%%%%%%
%         h11 = h11 + abs(y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,1))*xx_r(k,:)'*xx_r(k,:);
%         h22 = h22 + abs(y(k,2)-theta(k).*y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,2))*xx_r(k,:)'*xx_r(k,:);
%         h33 = h33 + exp(-xx_r(k,:)*gamma_1(:,2))*ck*xx_r(k,:)'*xx_r(k,:);
%         h44 = h44 + exp(-xx_r(k,:)*gamma_1(:,2))*ck*xx_i(k,:)'*xx_i(k,:);
%         h42 = h42 + exp(-xx_r(k,:)*gamma_1(:,2))*(ck*(xx_i(k,:)*gamma_2(:,1))*xx_i(k,:)+ik*xx_i(k,:))'*xx_r(k,:);
%         h32 = h32 + exp(-xx_r(k,:)*gamma_1(:,2))*(ck*(xx_r(k,:)*gamma_1(:,3))*xx_r(k,:)+rk*xx_r(k,:))'*xx_r(k,:);
%     end
%     
%     %first frequency
%     rk = -y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1));
%     ik = sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1)));
%     ck = 2*abs(y(1,1)).^2;
%     gr1 = gr1 + 0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))));
%     gr2 = gr2 + 0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))));  
%     gr3 = gr3 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*rk*xx_r(1,:)' +...
%                 exp(-xx_r(1,:)*gamma_1(:,2))*ck*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)');
%     gr4 = gr4 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*ik*xx_i(1,:)' +...
%                 exp(-xx_r(1,:)*gamma_1(:,2))*ck*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)');
%              
%     h11 = h11 + 0.5*(abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))*xx_r(1,:)'*xx_r(1,:));
%     h22 = h22 + 0.5*(abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))*xx_r(1,:)'*xx_r(1,:));
%     h33 = h33 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*ck*xx_r(1,:)'*xx_r(1,:));
%     h44 = h44 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*ck*xx_i(1,:)'*xx_i(1,:));
%     h42 = h42 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(ck*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)+ik*xx_i(1,:))'*xx_r(1,:));
%     h32 = h32 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(ck*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)+rk*xx_r(1,:))'*xx_r(1,:));
    h24=h42';
    h23=h32';
    
else
    %%%%%%%%%%%%%%%%%%%%%
    %gradient
    %%%%%%%%%%%%%%%%%%%%%
    rk = -y(2:nfreq,1).*conj(y(2:nfreq,2)) - y(2:nfreq,2).*conj(y(2:nfreq,1));
    ik = sqrt(-1)*(-y(2:nfreq,1).*conj(y(2:nfreq,2)) + y(2:nfreq,2).*conj(y(2:nfreq,1)));
    ck = 2*abs(y(2:nfreq,1)).^2;
    
    gr1 = gr1 + xx_r(2:nfreq,:)'*(1-abs(y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*gamma_1(:,1))) + ...
              0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1)))) +...
              0.5*(xx_r(end,:)'*(1-abs(y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,1))));
    gr2 = gr2 + xx_r(2:nfreq,:)'*(1 - abs(y(2:nfreq,2)-theta(2:nfreq).*y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*gamma_1(:,2))) + ...
              0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2)))) + ...
              0.5*(xx_r(end,:)'*(1 - abs(y(end,2)-theta(end).*y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,2))));  
    temp_mat_31 = bsxfun(@times,rk, xx_r(2:nfreq,:));
    temp_mat_32 = bsxfun(@times,ck,bsxfun(@times, xx_r(2:nfreq,:),xx_r(2:nfreq,:)*gamma_1(:,3)));
    gr3 = gr3 + sum( bsxfun(@times, (temp_mat_31 + temp_mat_32), exp(-xx_r(2:nfreq,:)*gamma_1(:,2)) ))' +...
          0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:)' +...
          exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)') +...
          0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*(-y(end,1).*conj(y(end,2)) - y(end,2).*conj(y(end,1)))*xx_r(end,:)' +...
          exp(-xx_r(end,:)*gamma_1(:,2))*2*abs(y(end,1)).^2*(xx_r(end,:)*gamma_1(:,3))*xx_r(end,:)');
      
    temp_mat_41 = bsxfun(@times, ik, xx_i(2:nfreq,:));
    temp_mat_42 = bsxfun(@times, ck, bsxfun(@times,xx_i(2:nfreq,:),xx_i(2:nfreq,:)*gamma_2(:,1)));
    gr4 = gr4 + sum( bsxfun(@times, (temp_mat_41 + temp_mat_42), exp(-xx_r(2:nfreq,:)*gamma_1(:,2))))' + ...
          0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1))))*xx_i(1,:)' +...
          exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)') +...
          0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*(sqrt(-1)*(-y(end,1).*conj(y(end,2)) + y(end,2).*conj(y(end,1))))*xx_i(end,:)' +...
          exp(-xx_r(end,:)*gamma_1(:,2))*2*abs(y(end,1)).^2*(xx_i(end,:)*gamma_2(:,1))*xx_i(end,:)');
    %%%%%%%%%%%%%%%%%%%
    %hession
    %%%%%%%%%%%%%%%%%%%
    bigmat_h11 = kron(bsxfun(@times, abs(y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*gamma_1(:,1)), xx_r(2:nfreq,:)),...
                                                                                ones(ngamma,1)');
    coefmat_h11 = repmat(xx_r(2:nfreq,:), 1,ngamma);
    h11 = h11 + reshape(sum(bsxfun(@times, bigmat_h11, coefmat_h11),1),ngamma,ngamma) +...
            0.5*(abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))*xx_r(1,:)'*xx_r(1,:))+...
            0.5*(abs(y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,1))*xx_r(end,:)'*xx_r(end,:));
    
    bigmat_h22 = kron(bsxfun(@times,  abs(y(2:nfreq,2)-theta(2:nfreq).*y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*gamma_1(:,2)),...
                                                                xx_r(2:nfreq,:)), ones(ngamma,1)');
    coefmat_h22 = repmat(xx_r(2:nfreq,:), 1,ngamma);
    h22 = h22 + reshape(sum(bsxfun(@times, bigmat_h22, coefmat_h22),1),ngamma,ngamma) +...
            0.5*(abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))*xx_r(1,:)'*xx_r(1,:))+...
            0.5*(abs(y(end,2)-theta(end).*y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,2))*xx_r(end,:)'*xx_r(end,:));
    
    bigmat_h33 = kron(bsxfun(@times, exp(-xx_r(2:nfreq,:)*gamma_1(:,2)).*ck, xx_r(2:nfreq,:)),...
                                                                                ones(ngamma,1)');
    coefmat_h33 = repmat(xx_r(2:nfreq,:), 1,ngamma);         
    h33 = h33 + reshape(sum(bsxfun(@times, bigmat_h33, coefmat_h33),1),ngamma,ngamma) +...
                0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*xx_r(1,:)'*xx_r(1,:))+...
                0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*2*abs(y(end,1)).^2*xx_r(end,:)'*xx_r(end,:));


    bigmat_h44 = kron(bsxfun(@times, exp(-xx_r(2:nfreq,:)*gamma_1(:,2)).*ck, xx_i(2:nfreq,:)),...
                                                                                ones(nbasis,1)');
    coefmat_h44 = repmat(xx_i(2:nfreq,:), 1,nbasis);         
    h44 = h44 + reshape(sum(bsxfun(@times, bigmat_h44, coefmat_h44),1),nbasis,nbasis) +...
                0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*2*abs(y(1,1)).^2*xx_i(1,:)'*xx_i(1,:))+...
                0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*2*abs(y(end,1)).^2*xx_i(end,:)'*xx_i(end,:));
    
    bigmat_h42_1 = kron(bsxfun(@times, exp(-xx_r(2:nfreq,:)*gamma_1(:,2)).*ck.*(xx_i(2:nfreq,:)*gamma_2(:,1)),...
                        xx_i(2:nfreq,:)), ones(ngamma,1)');
    coefmat_h42_1 = repmat(xx_r(2:nfreq,:), 1,nbasis);  
    bigmat_h42_2 = kron(bsxfun(@times,exp(-xx_r(2:nfreq,:)*gamma_1(:,2)).*ik, xx_i(2:nfreq,:)), ones(ngamma,1)');
    coefmat_h42_2 = repmat(xx_r(2:nfreq,:), 1,nbasis);  
    h42 = h42 + reshape(sum(bsxfun(@times, bigmat_h42_1, coefmat_h42_1) +...
                        bsxfun(@times, bigmat_h42_2, coefmat_h42_2),1),ngamma,nbasis)' +...
                        0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(2*abs(y(1,1)).^2*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)+...
                        sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1)))*xx_i(1,:))'*xx_r(1,:)) +...
                        0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*(2*abs(y(end,1)).^2*(xx_i(end,:)*gamma_2(:,1))*xx_i(end,:)+...
                        sqrt(-1)*(-y(end,1).*conj(y(end,2)) + y(end,2).*conj(y(end,1)))*xx_i(end,:))'*xx_r(end,:));

    bigmat_h32_1 = kron(bsxfun(@times, exp(-xx_r(2:nfreq,:)*gamma_1(:,2)).*ck.*(xx_r(2:nfreq,:)*gamma_1(:,3)),...
                        xx_r(2:nfreq,:)), ones(ngamma,1)');
    coefmat_h32_1 = repmat(xx_r(2:nfreq,:), 1,ngamma);  
    bigmat_h32_2 = kron(bsxfun(@times,exp(-xx_r(2:nfreq,:)*gamma_1(:,2)).*rk, xx_r(2:nfreq,:)), ones(ngamma,1)');
    coefmat_h32_2 = repmat(xx_r(2:nfreq,:), 1,ngamma);  
    h32 = h32 + reshape(sum(bsxfun(@times, bigmat_h32_1, coefmat_h32_1) +...
                        bsxfun(@times, bigmat_h32_2, coefmat_h32_2),1),ngamma,ngamma)' +...
                        0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(2*abs(y(1,1)).^2*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)+...
                        (-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:))'*xx_r(1,:)) +...
                        0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*(2*abs(y(end,1)).^2*(xx_r(end,:)*gamma_1(:,3))*xx_r(end,:)+...
                        (-y(end,1).*conj(y(end,2)) - y(end,2).*conj(y(end,1)))*xx_r(end,:))'*xx_r(end,:));
 
%     for k=2:nfreq
%         
%         rk = -y(k,1).*conj(y(k,2)) - y(k,2).*conj(y(k,1));
%         ik = sqrt(-1)*(-y(k,1).*conj(y(k,2)) + y(k,2).*conj(y(k,1)));
%         ck = 2*abs(y(k,1)).^2;
%         gr1 = gr1 + xx_r(k,:)'*(1-abs(y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,1)));
%         gr2 = gr2 + xx_r(k,:)'*(1 - abs(y(k,2)-theta(k).*y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,2)));  
%         gr3 = gr3 + exp(-xx_r(k,:)*gamma_1(:,2))*rk*xx_r(k,:)' +...
%                     exp(-xx_r(k,:)*gamma_1(:,2))*ck*(xx_r(k,:)*gamma_1(:,3))*xx_r(k,:)';
%         gr4 = gr4 + exp(-xx_r(k,:)*gamma_1(:,2))*ik*xx_i(k,:)' +...
%                     exp(-xx_r(k,:)*gamma_1(:,2))*ck*(xx_i(k,:)*gamma_2(:,1))*xx_i(k,:)';
% 
%         %%%%%%%%%%%%%%%%%%%%% 
%         %hession        
%         %%%%%%%%%%%%%%%%%%%%%
%         h11 = h11 + abs(y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,1))*xx_r(k,:)'*xx_r(k,:);
%         h22 = h22 + abs(y(k,2)-theta(k).*y(k,1)).^2.*exp(-xx_r(k,:)*gamma_1(:,2))*xx_r(k,:)'*xx_r(k,:);
%         h33 = h33 + exp(-xx_r(k,:)*gamma_1(:,2))*ck*xx_r(k,:)'*xx_r(k,:);
%         h44 = h44 + exp(-xx_r(k,:)*gamma_1(:,2))*ck*xx_i(k,:)'*xx_i(k,:);
%         h42 = h42 + exp(-xx_r(k,:)*gamma_1(:,2))*(ck*(xx_i(k,:)*gamma_2(:,1))*xx_i(k,:)+ik*xx_i(k,:))'*xx_r(k,:);
%         h32 = h32 + exp(-xx_r(k,:)*gamma_1(:,2))*(ck*(xx_r(k,:)*gamma_1(:,3))*xx_r(k,:)+rk*xx_r(k,:))'*xx_r(k,:);
%     end 
%     
%     %first frequency
%     rk = -y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1));
%     ik = sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1)));
%     ck = 2*abs(y(1,1)).^2;
%     gr1 = gr1 + 0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))));
%     gr2 = gr2 + 0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))));  
%     gr3 = gr3 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*rk*xx_r(1,:)' +...
%                 exp(-xx_r(1,:)*gamma_1(:,2))*ck*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)');
%     gr4 = gr4 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*ik*xx_i(1,:)' +...
%                 exp(-xx_r(1,:)*gamma_1(:,2))*ck*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)');
%               
%     h11 = h11 + 0.5*(abs(y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,1))*xx_r(1,:)'*xx_r(1,:));
%     h22 = h22 + 0.5*(abs(y(1,2)-theta(1).*y(1,1)).^2.*exp(-xx_r(1,:)*gamma_1(:,2))*xx_r(1,:)'*xx_r(1,:));
%     h33 = h33 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*ck*xx_r(1,:)'*xx_r(1,:));
%     h44 = h44 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*ck*xx_i(1,:)'*xx_i(1,:));
%     h42 = h42 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(ck*(xx_i(1,:)*gamma_2(:,1))*xx_i(1,:)+ik*xx_i(1,:))'*xx_r(1,:));
%     h32 = h32 + 0.5*(exp(-xx_r(1,:)*gamma_1(:,2))*(ck*(xx_r(1,:)*gamma_1(:,3))*xx_r(1,:)+rk*xx_r(1,:))'*xx_r(1,:));
% 
%     %end frequency
%     rk = -y(end,1).*conj(y(end,2)) - y(end,2).*conj(y(end,1));
%     ik = sqrt(-1)*(-y(end,1).*conj(y(end,2)) + y(end,2).*conj(y(end,1)));
%     ck = 2*abs(y(end,1)).^2;
%     gr1 = gr1 + 0.5*(xx_r(end,:)'*(1-abs(y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,1))));
%     gr2 = gr2 + 0.5*(xx_r(end,:)'*(1 - abs(y(end,2)-theta(end).*y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,2))));  
%     gr3 = gr3 + 0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*rk*xx_r(end,:)' +...
%                 exp(-xx_r(end,:)*gamma_1(:,2))*ck*(xx_r(end,:)*gamma_1(:,3))*xx_r(end,:)');
%     gr4 = gr4 + 0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*ik*xx_i(end,:)' +...
%                 exp(-xx_r(end,:)*gamma_1(:,2))*ck*(xx_i(end,:)*gamma_2(:,1))*xx_i(end,:)');
%              
%     h11 = h11 + 0.5*(abs(y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,1))*xx_r(end,:)'*xx_r(end,:));
%     h22 = h22 + 0.5*(abs(y(end,2)-theta(end).*y(end,1)).^2.*exp(-xx_r(end,:)*gamma_1(:,2))*xx_r(end,:)'*xx_r(end,:));
%     h33 = h33 + 0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*ck*xx_r(end,:)'*xx_r(end,:));
%     h44 = h44 + 0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*ck*xx_i(end,:)'*xx_i(end,:));
%     h42 = h42 + 0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*(ck*(xx_i(end,:)*gamma_2(:,1))*xx_i(end,:)+ik*xx_i(end,:))'*xx_r(end,:));
%     h32 = h32 + 0.5*(exp(-xx_r(end,:)*gamma_1(:,2))*(ck*(xx_r(end,:)*gamma_1(:,3))*xx_r(end,:)+rk*xx_r(end,:))'*xx_r(end,:));
% 
    h24=h42';
    h23=h32';
end   

h1 = [h11,zeros(ngamma,2*ngamma+nbasis)];
h2 = [zeros(ngamma,ngamma),h22,-h23,-h24];
h3 = [zeros(ngamma,ngamma),-h32,h33,zeros(ngamma,nbasis)];
h4 = [zeros(nbasis,ngamma),-h42,zeros(nbasis,ngamma),h44];    

gr = [gr1;gr2;gr3;gr4];
h = [h1;h2;h3;h4];
f = -f; 

gr_index = (1:(3*ngamma+nbasis)).*[kron(chol_index(rho_temp,1:3),ones(ngamma,1)'), kron(chol_index(rho_temp,4),ones(nbasis,1)')];
gr_index = gr_index(find(gr_index~=0));
gr = gr(gr_index);
h = h(gr_index,gr_index);