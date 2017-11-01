function[spectra, coh] = MultiSpect_interval(zt, spect_matrices, fit, params, alphalevel) 

dim = size(zt); nobs = dim(1); dimen=dim(2);
spectra = cell(dimen,1);
coh = cell(dimen,1);
if dimen==3  
    
    temp1 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp2 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp3 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp4 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp5 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp6 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    spect11 = zeros(params.nfreq+1,nobs,3);
    spect22 = zeros(params.nfreq+1,nobs,3);
    spect33 = zeros(params.nfreq+1,nobs,3);
    coh21 = zeros(params.nfreq+1,nobs,3);
    coh31 = zeros(params.nfreq+1,nobs,3);
    coh32 = zeros(params.nfreq+1,nobs,3);
    for j=1:nobs
        if(mod(j,10)==0)
            fprintf('Completed: %g of %g observations \n' ,j, nobs)
        end
        for p=1:params.nloop
            if(p>params.nwarmup)
                if length(fit(fit(1).nexp_curr(p)).xi(:,p))==1
                    k=1;
                else    
                    k = min(find(j<=fit(fit(1).nexp_curr(p)).xi(:,p)));
                end    
                temp1(:,p-params.nwarmup) = squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,k,p));
                temp2(:,p-params.nwarmup) = squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,k,p));
                temp3(:,p-params.nwarmup) = squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,3,:,k,p));
                temp4(:,p-params.nwarmup) = abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,1,:,k,p))).^2./...
                                (squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,k,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,k,p)));
                temp5(:,p-params.nwarmup) = abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,1,:,k,p))).^2./...
                                (squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,k,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,3,:,k,p)));        
                temp6(:,p-params.nwarmup) = abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,2,:,k,p))).^2./...
                                (squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,k,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,3,:,k,p)));                         
            end
        end
        spect11(:,j,1) = mean(temp1,2); spect11(:,j,2:3) =  quantile(temp1,[alphalevel/2, 1-alphalevel/2],2);
        spect22(:,j,1) = mean(temp2,2); spect22(:,j,2:3) =  quantile(temp2,[alphalevel/2, 1-alphalevel/2],2);
        spect33(:,j,1) = mean(temp3,2); spect33(:,j,2:3) =  quantile(temp3,[alphalevel/2, 1-alphalevel/2],2);
        coh21(:,j,1) = mean(temp4,2); coh21(:,j,2:3) =  quantile(temp4,[alphalevel/2, 1-alphalevel/2],2);
        coh31(:,j,1) = mean(temp5,2); coh31(:,j,2:3) =  quantile(temp5,[alphalevel/2, 1-alphalevel/2],2);
        coh32(:,j,1) = mean(temp6,2); coh32(:,j,2:3) =  quantile(temp6,[alphalevel/2, 1-alphalevel/2],2);
    end    
    
    spectra{1}=spect11;spectra{2}=spect22;spectra{3}=spect33;
    coh{1}=coh21;coh{2}=coh31; coh{3}=coh32;
    
elseif dimen==2
  
    temp1 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp2 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    temp3 = zeros(params.nfreq+1,(params.nloop-params.nwarmup));
    spect11 = zeros(params.nfreq+1,nobs,3);
    spect22 = zeros(params.nfreq+1,nobs,3);
    coh21 = zeros(params.nfreq+1,nobs,3);
    for j=1:nobs
        if(mod(j,10)==0)
            fprintf('Completed: %g of %g observations \n' ,j, nobs)
        end
        for p=1:params.nloop
            if(p>params.nwarmup)
                if length(fit(fit(1).nexp_curr(p)).xi(:,p))==1
                    k=1;
                else    
                    k = min(find(j<=fit(fit(1).nexp_curr(p)).xi(:,p)));
                end    
                temp1(:,p-params.nwarmup) = squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,k,p));
                temp2(:,p-params.nwarmup) = squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,k,p));
                temp3(:,p-params.nwarmup) = abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,1,:,k,p))).^2./...
                                    (squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,k,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,k,p))); 
            end
        end
        spect11(:,j,1) = mean(temp1,2); spect11(:,j,2:3) =  quantile(temp1,[alphalevel/2, 1-alphalevel/2],2);
        spect22(:,j,1) = mean(temp2,2); spect22(:,j,2:3) =  quantile(temp2,[alphalevel/2, 1-alphalevel/2],2);
        coh21(:,j,1) = mean(temp3,2); coh21(:,j,2:3) =  quantile(temp3,[alphalevel/2, 1-alphalevel/2],2);
    end    

    spectra{1}=spect11;spectra{2}=spect22;
    coh{1}=coh21;  

end









