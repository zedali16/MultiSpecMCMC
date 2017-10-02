function [gr] = Gradient1(yobs_tmp, chol_index, Phi_temp, tau_temp,...
                                    Beta_temp, sigmasqalpha, nbasis)
global dimen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate log gradients for coeeficients selected to be different
%
%   Input:
%       1) x - initial values for coeeficient of basis functions need to be
%       optimized 
%       2) yobs_tmp - time series data within the segment
%       3) chol_index - index matrix
%       4) Phi_temp - which component changed
%       5) tau_temp - smoothing parameters  
%       6) Beta_temp - current coefficients
%       7) sigmasqalpha - smoothing parameters for the constant in real
%       components
%       8) nbasis - number of basis function used
%   Main Outputs:
%       2) gr - gradients for optimization process
%
%   Required programs: lin_basis_func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initilize Beta_1 and Beta_2: Beta_2 is for imaginary components
nBeta = nbasis + 1;

Beta_1(:,:) = Beta_temp(:,1:(dimen + dimen*(dimen-1)/2));
Beta_2(:,:) = Beta_temp(1:nBeta,(dimen + dimen*(dimen-1)/2 + 1): end);

dim = size(yobs_tmp); n = dim(1);
nfreq = floor(n/2); tt = (0:nfreq)/(2*nfreq);
yy = fft(yobs_tmp)/sqrt(n); y = yy(1:(nfreq+1),:); nf = length(y);
[xx_r, xx_i]=lin_basis_func(tt);

%theta's
theta = zeros(dimen*(dimen-1)/2,nf);
for i=1:dimen*(dimen-1)/2
    theta_real = xx_r * Beta_1(:,i+dimen);
    theta_imag = xx_i * Beta_2(:,i);
    theta(i,:) = theta_real + sqrt(-1)*theta_imag;
end    
%delta's
delta_sq = zeros(dimen,nf);
for i=1:dimen
    delta_sq(i,:) = exp(xx_r * Beta_1(:,i));
end   

if dimen==2  %Bivariate Time Series

    gr1 = zeros(nBeta,1); gr2 = zeros(nBeta,1); gr3 = zeros(nBeta,1); gr4 = zeros(nBeta,1);
    gr1(1) = Beta_1(1,1)/sigmasqalpha; gr1(2:nBeta,1) = Beta_1(2:nBeta,1)/tau_temp(1);
    gr2(1) = Beta_1(1,2)/sigmasqalpha; gr2(2:nBeta,1) = Beta_1(2:nBeta,2)/tau_temp(2);
    gr3(1) = Beta_1(1,3)/sigmasqalpha; gr3(2:nBeta,1) = Beta_1(2:nBeta,3)/tau_temp(3);
    gr4(1:nBeta,1) = Beta_2(1:nBeta,1)/tau_temp(4);

    if (mod(n,2)==1)     
        %%%%%%%%%%%%%%%%%%%%%%%%
        %gradient
        %%%%%%%%%%%%%%%%%%%%%%%%
        rk = -y(2:end,1).*conj(y(2:end,2)) - y(2:end,2).*conj(y(2:end,1));
        ik = sqrt(-1)*(-y(2:end,1).*conj(y(2:end,2)) + y(2:end,2).*conj(y(2:end,1)));
        ck = 2*abs(y(2:end,1)).^2;
    
        gr1 = gr1 + xx_r(2:end,:)'*(1-abs(y(2:end,1)).^2.*exp(-xx_r(2:end,:)*Beta_1(:,1))) + ...
                    0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,1))));
        gr2 = gr2 + xx_r(2:end,:)'*(1 - abs(y(2:end,2)-theta(2:end).'.*y(2:end,1)).^2.*exp(-xx_r(2:end,:)*Beta_1(:,2))) + ...
                    0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1).'.*y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,2))));  
        temp_mat_31 = bsxfun(@times, xx_r(2:end,:),rk);
        temp_mat_32 = bsxfun(@times, ck, bsxfun(@times, xx_r(2:end,:),xx_r(2:end,:)*Beta_1(:,3)));
        gr3 = gr3 + sum( bsxfun(@times, (temp_mat_31+temp_mat_32), exp(-xx_r(2:end,:)*Beta_1(:,2)) ))' +...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_r(1,:)*Beta_1(:,3))*xx_r(1,:)');
        temp_mat_41 = bsxfun(@times, ik, xx_i(2:end,:));
        temp_mat_42 = bsxfun(@times, ck, bsxfun(@times,xx_i(2:end,:),xx_i(2:end,:)*Beta_2(:,1)));
        gr4 = gr4 + sum( bsxfun(@times, (temp_mat_41 + temp_mat_42), exp(-xx_r(2:end,:)*Beta_1(:,2))))' + ...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_i(1,:)*Beta_2(:,1))*xx_i(1,:)');    
    else
        %%%%%%%%%%%%%%%%%%%%%
        %gradient
        %%%%%%%%%%%%%%%%%%%%%
        rk = -y(2:nfreq,1).*conj(y(2:nfreq,2)) - y(2:nfreq,2).*conj(y(2:nfreq,1));
        ik = sqrt(-1)*(-y(2:nfreq,1).*conj(y(2:nfreq,2)) + y(2:nfreq,2).*conj(y(2:nfreq,1)));
        ck = 2*abs(y(2:nfreq,1)).^2;
    
        gr1 = gr1 + xx_r(2:nfreq,:)'*(1-abs(y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*Beta_1(:,1))) + ...
                    0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,1)))) +...
                    0.5*(xx_r(end,:)'*(1-abs(y(end,1)).^2.*exp(-xx_r(end,:)*Beta_1(:,1))));
        gr2 = gr2 + xx_r(2:nfreq,:)'*(1 - abs(y(2:nfreq,2)-theta(2:nfreq).'.*y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*Beta_1(:,2))) + ...
                    0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1).'.*y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,2)))) + ...
                    0.5*(xx_r(end,:)'*(1 - abs(y(end,2)-theta(end).'.*y(end,1)).^2.*exp(-xx_r(end,:)*Beta_1(:,2))));  
        temp_mat_31 = bsxfun(@times,rk, xx_r(2:nfreq,:));
        temp_mat_32 = bsxfun(@times,ck,bsxfun(@times, xx_r(2:nfreq,:),xx_r(2:nfreq,:)*Beta_1(:,3)));
        gr3 = gr3 + sum( bsxfun(@times, (temp_mat_31 + temp_mat_32), exp(-xx_r(2:nfreq,:)*Beta_1(:,2)) ))' +...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_r(1,:)*Beta_1(:,3))*xx_r(1,:)') +...
                    0.5*(exp(-xx_r(end,:)*Beta_1(:,2))*(-y(end,1).*conj(y(end,2)) - y(end,2).*conj(y(end,1)))*xx_r(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,2))*2*abs(y(end,1)).^2*(xx_r(end,:)*Beta_1(:,3))*xx_r(end,:)');
        temp_mat_41 = bsxfun(@times, ik, xx_i(2:nfreq,:));
        temp_mat_42 = bsxfun(@times, ck, bsxfun(@times,xx_i(2:nfreq,:),xx_i(2:nfreq,:)*Beta_2(:,1)));
        gr4 = gr4 + sum( bsxfun(@times, (temp_mat_41 + temp_mat_42), exp(-xx_r(2:nfreq,:)*Beta_1(:,2))))' + ...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(sqrt(-1)*(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_i(1,:)*Beta_2(:,1))*xx_i(1,:)') +...
                    0.5*(exp(-xx_r(end,:)*Beta_1(:,2))*(sqrt(-1)*(-y(end,1).*conj(y(end,2)) + y(end,2).*conj(y(end,1))))*xx_i(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,2))*2*abs(y(end,1)).^2*(xx_i(end,:)*Beta_2(:,1))*xx_i(end,:)');
    end   
    gr = [gr1;gr2;gr3;gr4];
    gr_index = (1:(4*nBeta)).*[kron(chol_index(Phi_temp,1:3),ones(nBeta,1)'), kron(chol_index(Phi_temp,4),ones(nBeta,1)')];
    gr_index = gr_index(find(gr_index~=0));
    gr = gr(gr_index); 
    
elseif dimen==3  %trivariate time series
    
    gr1 = zeros(nBeta,1); gr2 = zeros(nBeta,1); gr3 = zeros(nBeta,1); gr4 = zeros(nBeta,1);
    gr5 = zeros(nBeta,1); gr6 = zeros(nBeta,1); gr7 = zeros(nBeta,1); gr8 = zeros(nBeta,1);
    gr9 = zeros(nBeta,1);

    gr1(1) = Beta_1(1,1)/sigmasqalpha; gr1(2:nBeta) = Beta_1(2:nBeta,1)/tau_temp(1);
    gr2(1) = Beta_1(1,2)/sigmasqalpha; gr2(2:nBeta) = Beta_1(2:nBeta,2)/tau_temp(2);
    gr3(1) = Beta_1(1,3)/sigmasqalpha; gr3(2:nBeta) = Beta_1(2:nBeta,3)/tau_temp(3);
    gr4(1) = Beta_1(1,4)/sigmasqalpha; gr4(2:nBeta) = Beta_1(2:nBeta,4)/tau_temp(4);
    gr5(1) = Beta_1(1,5)/sigmasqalpha; gr5(2:nBeta) = Beta_1(2:nBeta,5)/tau_temp(5);
    gr6(1) = Beta_1(1,6)/sigmasqalpha; gr6(2:nBeta) = Beta_1(2:nBeta,6)/tau_temp(6);
    gr7(1:nBeta) = Beta_2(1:nBeta,1)/tau_temp(7);
    gr8(1:nBeta) = Beta_2(1:nBeta,2)/tau_temp(8);
    gr9(1:nBeta) = Beta_2(1:nBeta,3)/tau_temp(9);
    
    if (mod(n,2)==1)     
        %%%%%%%%%%%%%%%%%%%%%%%%
        %gradient
        %%%%%%%%%%%%%%%%%%%%%%%%
        rk4 = -y(2:end,1).*conj(y(2:end,2)) - y(2:end,2).*conj(y(2:end,1));
        ck4 = 2*abs(y(2:end,1)).^2;
        rk5 = -y(2:end,1).*conj(y(2:end,3)) - y(2:end,3).*conj(y(2:end,1));
        ck5 = 2*abs(y(2:end,1)).^2;
        b = theta(3,:);
        dk5 =  y(2:end,2).*conj(y(2:end,1)).*(b(2:end).') + conj(y(2:end,2)).*y(2:end,1).*conj(b(2:end).');
        rk6 = -y(2:end,2).*conj(y(2:end,3)) - y(2:end,3).*conj(y(2:end,2));
        ck6 = 2*abs(y(2:end,2)).^2;
        a = theta(2,:);
        dk6 =  y(2:end,1).*conj(y(2:end,2)).*(a(2:end).') + conj(y(2:end,1)).*y(2:end,2).*conj(a(2:end).');
        ik7 = sqrt(-1)*(-y(2:end,1).*conj(y(2:end,2)) + y(2:end,2).*conj(y(2:end,1)));
        ck7 = 2*abs(y(2:end,1)).^2;
        ik8 = sqrt(-1)*(-y(2:end,1).*conj(y(2:end,3)) + y(2:end,3).*conj(y(2:end,1)));
        ck8 = 2*abs(y(2:end,1)).^2;
        dk8 = sqrt(-1)*(-y(2:end,2).*conj(y(2:end,1)).*(b(2:end).') + conj(y(2:end,2).*conj(y(2:end,1)).*(b(2:end).'))) ;
        ik9 = sqrt(-1)*(-y(2:end,2).*conj(y(2:end,3)) + y(2:end,3).*conj(y(2:end,2)));
        ck9 = 2*abs(y(2:end,2)).^2;
        dk9 = sqrt(-1)*(-y(2:end,1).*conj(y(2:end,2)).*(a(2:end).') + conj(y(2:end,1).*conj(y(2:end,2)).*(a(2:end).'))) ;
    
        gr1 = gr1 + xx_r(2:end,:)'*(1-abs(y(2:end,1)).^2.*exp(-xx_r(2:end,:)*Beta_1(:,1))) + ...
                    0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,1))));
        gr2 = gr2 + xx_r(2:end,:)'*(1 - abs(y(2:end,2)-theta(1,2:end).'.*y(2:end,1)).^2.*exp(-xx_r(2:end,:)*Beta_1(:,2))) + ...
                    0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1,1).*y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,2))));  
        gr3 = gr3 + xx_r(2:end,:)'*(1 - abs(y(2:end,3) - theta(2,2:end).'.*y(2:end,1) - theta(3,2:end).'.*y(2:end,2)).^2.*exp(-xx_r(2:end,:)*Beta_1(:,3)))+...
                    0.5*(xx_r(1,:)'*(1 - abs(y(1,3) - theta(2,1).*y(1,1) - theta(3,1).*y(1,2))^2.*exp(-xx_r(1,:)*Beta_1(:,3))));             
        temp_mat_41 = bsxfun(@times, xx_r(2:end,:),rk4);
        temp_mat_42 = bsxfun(@times, ck4, bsxfun(@times, xx_r(2:end,:),xx_r(2:end,:)*Beta_1(:,4)));
        gr4 = gr4 + sum( bsxfun(@times, (temp_mat_41+temp_mat_42), exp(-xx_r(2:end,:)*Beta_1(:,2)) ))' +...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_r(1,:)*Beta_1(:,4))*xx_r(1,:)');
        temp_mat_51 = bsxfun(@times, xx_r(2:end,:),rk5);
        temp_mat_52 = bsxfun(@times, ck5, bsxfun(@times, xx_r(2:end,:),xx_r(2:end,:)*Beta_1(:,5))); 
        temp_mat_53 = bsxfun(@times, xx_r(2:end,:),dk5);
        gr5 = gr5 + sum( bsxfun(@times, (temp_mat_51 + temp_mat_52 + temp_mat_53), exp(-xx_r(2:end,:)*Beta_1(:,3)) ))'+...
                    0.5* ( exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,1).*conj(y(1,3)) - y(1,3).*conj(y(1,1)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,1)).^2*(xx_r(1,:)*Beta_1(:,5))*xx_r(1,:)'+...
                    exp(-xx_r(1,:)*Beta_1(:,3))*(y(1,2).*conj(y(1,1)).*(b(1).') + conj(y(1,2).*conj(y(1,1)).*(b(1).')))*xx_r(1,:)');
        temp_mat_61 = bsxfun(@times, xx_r(2:end,:),rk6);
        temp_mat_62 = bsxfun(@times, ck6, bsxfun(@times, xx_r(2:end,:),xx_r(2:end,:)*Beta_1(:,6))); 
        temp_mat_63 = bsxfun(@times, xx_r(2:end,:),dk6);
        gr6 = gr6 + sum(bsxfun(@times, (temp_mat_61 + temp_mat_62 + temp_mat_63), exp(-xx_r(2:end,:)*Beta_1(:,3)) ))'+...
                    0.5* (exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,2).*conj(y(1,3)) - y(1,3).*conj(y(1,2)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,2)).^2*(xx_r(1,:)*Beta_1(:,6))*xx_r(1,:)'+...
                    exp(-xx_r(1,:)*Beta_1(:,3))*(y(1,1).*conj(y(1,2)).*(a(1).') + conj(y(1,1).*conj(y(1,2)).*(a(1).')))*xx_r(1,:)');
        temp_mat_71 = bsxfun(@times, ik7, xx_i(2:end,:));
        temp_mat_72 = bsxfun(@times, ck7, bsxfun(@times,xx_i(2:end,:),xx_i(2:end,:)*Beta_2(:,1)));
        gr7 = gr7 + sum( bsxfun(@times, (temp_mat_71 + temp_mat_72), exp(-xx_r(2:end,:)*Beta_1(:,2))))' + ...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(imag(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_i(1,:)*Beta_2(:,1))*xx_i(1,:)');
        temp_mat_81 = bsxfun(@times, ik8, xx_i(2:end,:));
        temp_mat_82 = bsxfun(@times, ck8, bsxfun(@times,xx_i(2:end,:),xx_i(2:end,:)*Beta_2(:,2)));      
        temp_mat_83 = bsxfun(@times, xx_i(2:end,:),dk8);
        gr8 = gr8 + sum( bsxfun(@times, (temp_mat_81 + temp_mat_82 + temp_mat_83), exp(-xx_r(2:end,:)*Beta_1(:,3))))' + ...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,3))*(imag(-y(1,1).*conj(y(1,3)) + y(1,3).*conj(y(1,1))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,1)).^2*(xx_i(1,:)*Beta_2(:,2))*xx_i(1,:)'+...
                    sqrt(-1)*exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,1).*conj(y(1,2)).*(b(1).') + conj(y(1,1).*conj(y(1,2)).*(b(1).')))*xx_i(1,:)');
        temp_mat_91 = bsxfun(@times, ik9, xx_i(2:end,:));
        temp_mat_92 = bsxfun(@times, ck9, bsxfun(@times,xx_i(2:end,:),xx_i(2:end,:)*Beta_2(:,3)));      
        temp_mat_93 = bsxfun(@times, xx_i(2:end,:),dk9);
        gr9 = gr9 + sum( bsxfun(@times, (temp_mat_91 + temp_mat_92 + temp_mat_93), exp(-xx_r(2:end,:)*Beta_1(:,3))))' + ...
                    0.5 * (exp(-xx_r(1,:)*Beta_1(:,3))*(imag(-y(1,2).*conj(y(1,3)) + y(1,3).*conj(y(1,2))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,2)).^2*(xx_i(1,:)*Beta_2(:,3))*xx_i(1,:)'+...
                    sqrt(-1)*exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,1).*conj(y(1,2)).*(a(1).') + conj(y(1,1).*conj(y(1,2)).*(a(1).')))*xx_i(1,:)');     
    else
        %%%%%%%%%%%%%%%%%%%%%%%%
        %gradient
        %%%%%%%%%%%%%%%%%%%%%%%%
        rk4 = -y(2:nfreq,1).*conj(y(2:nfreq,2)) - y(2:nfreq,2).*conj(y(2:nfreq,1));
        ck4 = 2*abs(y(2:nfreq,1)).^2;
        rk5 = -y(2:nfreq,1).*conj(y(2:nfreq,3)) - y(2:nfreq,3).*conj(y(2:nfreq,1));
        ck5 = 2*abs(y(2:nfreq,1)).^2;
        b = theta(3,:);
        dk5 =  y(2:nfreq,2).*conj(y(2:nfreq,1)).*(b(2:nfreq).') + conj(y(2:nfreq,2)).*y(2:nfreq,1).*conj(b(2:nfreq).');
        rk6 = -y(2:nfreq,2).*conj(y(2:nfreq,3)) - y(2:nfreq,3).*conj(y(2:nfreq,2));
        ck6 = 2*abs(y(2:nfreq,2)).^2;
        a = theta(2,:);
        dk6 =  y(2:nfreq,1).*conj(y(2:nfreq,2)).*(a(2:nfreq).') + conj(y(2:nfreq,1).*conj(y(2:nfreq,2)).*(a(2:nfreq).'));
        ik7 = sqrt(-1)*(-y(2:nfreq,1).*conj(y(2:nfreq,2)) + y(2:nfreq,2).*conj(y(2:nfreq,1)));
        ck7 = 2*abs(y(2:nfreq,1)).^2;
        ik8 = sqrt(-1)*(-y(2:nfreq,1).*conj(y(2:nfreq,3)) + y(2:nfreq,3).*conj(y(2:nfreq,1)));
        ck8 = 2*abs(y(2:nfreq,1)).^2;
        dk8 = sqrt(-1)*( -y(2:nfreq,2).*conj(y(2:nfreq,1)).*(b(2:nfreq).') + conj(y(2:nfreq,2).*conj(y(2:nfreq,1)).*(b(2:nfreq).')));
        ik9 = sqrt(-1)*(-y(2:nfreq,2).*conj(y(2:nfreq,3)) + y(2:nfreq,3).*conj(y(2:nfreq,2)));
        ck9 = 2*abs(y(2:nfreq,2)).^2;
        dk9 = sqrt(-1)*( -y(2:nfreq,1).*conj(y(2:nfreq,2)).*(a(2:nfreq).') + conj(y(2:nfreq,1).*conj(y(2:nfreq,2)).*(a(2:nfreq).')));
    
        gr1 = gr1 + xx_r(2:nfreq,:)'*(1-abs(y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*Beta_1(:,1))) + ...
                    0.5*(xx_r(1,:)'*(1-abs(y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,1))))+...
                    0.5*(xx_r(end,:)'*(1-abs(y(end,1)).^2.*exp(-xx_r(end,:)*Beta_1(:,1))));
        gr2 = gr2 + xx_r(2:nfreq,:)'*(1 - abs(y(2:nfreq,2)-theta(1,2:nfreq).'.*y(2:nfreq,1)).^2.*exp(-xx_r(2:nfreq,:)*Beta_1(:,2))) + ...
                    0.5*(xx_r(1,:)'*(1 - abs(y(1,2)-theta(1,1).*y(1,1)).^2.*exp(-xx_r(1,:)*Beta_1(:,2)))) +...
                    0.5*(xx_r(end,:)'*(1 - abs(y(end,2)-theta(1,end).*y(end,1)).^2.*exp(-xx_r(end,:)*Beta_1(:,2))));  
        gr3 = gr3 + xx_r(2:nfreq,:)'*(1 - abs(y(2:nfreq,3) - theta(2,2:nfreq).'.*y(2:nfreq,1) - theta(3,2:nfreq).'.*y(2:nfreq,2)).^2.*exp(-xx_r(2:nfreq,:)*Beta_1(:,3)))+...
                    0.5*(xx_r(1,:)'*(1 - abs(y(1,3) - theta(2,1).*y(1,1) - theta(3,1).*y(1,2))^2.*exp(-xx_r(1,:)*Beta_1(:,3))))+...
                    0.5*(xx_r(end,:)'*(1 - abs(y(end,3) - theta(2,end).*y(end,1) - theta(3,end).*y(end,2))^2.*exp(-xx_r(end,:)*Beta_1(:,3))));             
        temp_mat_41 = bsxfun(@times, xx_r(2:nfreq,:),rk4);
        temp_mat_42 = bsxfun(@times, ck4, bsxfun(@times, xx_r(2:nfreq,:),xx_r(2:nfreq,:)*Beta_1(:,4)));
        gr4 = gr4 + sum( bsxfun(@times, (temp_mat_41+temp_mat_42), exp(-xx_r(2:nfreq,:)*Beta_1(:,2)) ))' +...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(-y(1,1).*conj(y(1,2)) - y(1,2).*conj(y(1,1)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_r(1,:)*Beta_1(:,4))*xx_r(1,:)')+...
                    0.5*(exp(-xx_r(end,:)*Beta_1(:,2))*(-y(end,1).*conj(y(end,2)) - y(end,2).*conj(y(end,1)))*xx_r(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,2))*2*abs(y(end,1)).^2*(xx_r(end,:)*Beta_1(:,4))*xx_r(end,:)');
        temp_mat_51 = bsxfun(@times, xx_r(2:nfreq,:),rk5);
        temp_mat_52 = bsxfun(@times, ck5, bsxfun(@times, xx_r(2:nfreq,:),xx_r(2:nfreq,:)*Beta_1(:,5))); 
        temp_mat_53 = bsxfun(@times, xx_r(2:nfreq,:),dk5);
        gr5 = gr5 + sum( bsxfun(@times, (temp_mat_51 + temp_mat_52 + temp_mat_53), exp(-xx_r(2:nfreq,:)*Beta_1(:,3)) ))'+...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,1).*conj(y(1,3)) - y(1,3).*conj(y(1,1)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,1)).^2*(xx_r(1,:)*Beta_1(:,5))*xx_r(1,:)'+...
                    exp(-xx_r(1,:)*Beta_1(:,3))*(y(1,2).*conj(y(1,1)).*(b(1).') + conj(y(1,2).*conj(y(1,1)).*(b(1).')))*xx_r(1,:)')+...
                    0.5*(exp(-xx_r(end,:)*Beta_1(:,3))*(-y(end,1).*conj(y(end,3)) - y(end,3).*conj(y(end,1)))*xx_r(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,3))*2*abs(y(end,1)).^2*(xx_r(end,:)*Beta_1(:,5))*xx_r(end,:)'+...
                    exp(-xx_r(end,:)*Beta_1(:,3))*(y(end,2).*conj(y(end,1)).*(b(end).') + conj(y(end,2).*conj(y(end,1)).*(b(end).')))*xx_r(end,:)');
        temp_mat_61 = bsxfun(@times, xx_r(2:nfreq,:),rk6);
        temp_mat_62 = bsxfun(@times, ck6, bsxfun(@times, xx_r(2:nfreq,:),xx_r(2:nfreq,:)*Beta_1(:,6))); 
        temp_mat_63 = bsxfun(@times, xx_r(2:nfreq,:),dk6);
        gr6 = gr6 + sum( bsxfun(@times, (temp_mat_61 + temp_mat_62 + temp_mat_63), exp(-xx_r(2:nfreq,:)*Beta_1(:,3)) ))'+...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,2).*conj(y(1,3)) - y(1,3).*conj(y(1,2)))*xx_r(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,2)).^2*(xx_r(1,:)*Beta_1(:,6))*xx_r(1,:)'+...
                    exp(-xx_r(1,:)*Beta_1(:,3))*(y(1,1).*conj(y(1,2)).*(a(1).') + conj(y(1,1).*conj(y(1,2)).*(a(1).')))*xx_r(1,:)')+...
                    0.5*(exp(-xx_r(end,:)*Beta_1(:,3))*(-y(end,2).*conj(y(end,3)) - y(end,3).*conj(y(end,2)))*xx_r(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,3))*2*abs(y(end,2)).^2*(xx_r(end,:)*Beta_1(:,6))*xx_r(end,:)'+...
                    exp(-xx_r(end,:)*Beta_1(:,3))*(y(end,1).*conj(y(end,2)).*(a(end).') + conj(y(end,1).*conj(y(end,2)).*(a(end).')))*xx_r(end,:)');
        temp_mat_71 = bsxfun(@times, ik7, xx_i(2:nfreq,:));
        temp_mat_72 = bsxfun(@times, ck7, bsxfun(@times,xx_i(2:nfreq,:),xx_i(2:nfreq,:)*Beta_2(:,1)));
        gr7 = gr7 + sum( bsxfun(@times, (temp_mat_71 + temp_mat_72), exp(-xx_r(2:nfreq,:)*Beta_1(:,2))))' + ...
                    0.5*(exp(-xx_r(1,:)*Beta_1(:,2))*(imag(-y(1,1).*conj(y(1,2)) + y(1,2).*conj(y(1,1))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,2))*2*abs(y(1,1)).^2*(xx_i(1,:)*Beta_2(:,1))*xx_i(1,:)')+...
                    0.5*(exp(-xx_r(end,:)*Beta_1(:,2))*(imag(-y(end,1).*conj(y(end,2)) + y(end,2).*conj(y(end,1))))*xx_i(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,2))*2*abs(y(end,1)).^2*(xx_i(end,:)*Beta_2(:,1))*xx_i(end,:)');
        temp_mat_81 = bsxfun(@times, ik8, xx_i(2:nfreq,:));
        temp_mat_82 = bsxfun(@times, ck8, bsxfun(@times,xx_i(2:nfreq,:),xx_i(2:nfreq,:)*Beta_2(:,2)));      
        temp_mat_83 = bsxfun(@times, xx_i(2:nfreq,:),dk8);
        gr8 = gr8 + sum( bsxfun(@times, (temp_mat_81 + temp_mat_82 + temp_mat_83), exp(-xx_r(2:nfreq,:)*Beta_1(:,3))))' + ...
                    0.5 * (exp(-xx_r(1,:)*Beta_1(:,3))*(imag(-y(1,1).*conj(y(1,3)) + y(1,3).*conj(y(1,1))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,1)).^2*(xx_i(1,:)*Beta_2(:,2))*xx_i(1,:)'+...
                    sqrt(-1)*exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,1).*conj(y(1,2)).*(b(1).') + conj(y(1,1).*conj(y(1,2)).*(b(1).')))*xx_i(1,:)')+...
                    0.5 * (exp(-xx_r(end,:)*Beta_1(:,3))*(imag(-y(end,1).*conj(y(end,3)) + y(end,3).*conj(y(end,1))))*xx_i(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,3))*2*abs(y(end,1)).^2*(xx_i(end,:)*Beta_2(:,2))*xx_i(end,:)'+...
                    sqrt(-1)*exp(-xx_r(end,:)*Beta_1(:,3))*(-y(end,1).*conj(y(end,2)).*(b(end).') + conj(y(end,1).*conj(y(end,2)).*(b(end).')))*xx_i(end,:)');
        temp_mat_91 = bsxfun(@times, ik9, xx_i(2:nfreq,:));
        temp_mat_92 = bsxfun(@times, ck9, bsxfun(@times,xx_i(2:nfreq,:),xx_i(2:nfreq,:)*Beta_2(:,3)));      
        temp_mat_93 = bsxfun(@times, xx_i(2:nfreq,:),dk9);
        gr9 = gr9 + sum( bsxfun(@times, (temp_mat_91 + temp_mat_92 + temp_mat_93), exp(-xx_r(2:nfreq,:)*Beta_1(:,3))))' + ...
                    0.5 * (exp(-xx_r(1,:)*Beta_1(:,3))*(imag(-y(1,2).*conj(y(1,3)) + y(1,3).*conj(y(1,2))))*xx_i(1,:)' +...
                    exp(-xx_r(1,:)*Beta_1(:,3))*2*abs(y(1,2)).^2*(xx_i(1,:)*Beta_2(:,3))*xx_i(1,:)'+...
                    sqrt(-1)*exp(-xx_r(1,:)*Beta_1(:,3))*(-y(1,1).*conj(y(1,2)).*(a(1).') + conj(y(1,1).*conj(y(1,2)).*(a(1).')))*xx_i(1,:)')+...
                    0.5 * (exp(-xx_r(end,:)*Beta_1(:,3))*(imag(-y(end,2).*conj(y(end,3)) + y(end,3).*conj(y(end,2))))*xx_i(end,:)' +...
                    exp(-xx_r(end,:)*Beta_1(:,3))*2*abs(y(end,2)).^2*(xx_i(end,:)*Beta_2(:,3))*xx_i(end,:)'+...
                    sqrt(-1)*exp(-xx_r(end,:)*Beta_1(:,3))*(-y(end,1).*conj(y(end,2)).*(a(end).') + conj(y(end,1).*conj(y(end,2)).*(a(end).')))*xx_i(end,:)');
    end
    gr = [gr1;gr2;gr3;gr4;gr5;gr6;gr7;gr8;gr9];
    gr_index = (1:(dimen^2*nBeta)).*kron(chol_index(Phi_temp,:),ones(nBeta,1)');
    gr_index = gr_index(find(gr_index~=0));
    gr = gr(gr_index); 
end    

gr=-gr;