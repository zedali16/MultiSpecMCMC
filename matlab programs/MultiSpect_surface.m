function[spectra, coh] = MultiSpect_surface(zt, spect_matrices, fit, params) 


dim = size(zt); nobs = dim(1); dimen=dim(2);
spectra = cell(dimen,1);
coh = cell(dimen,1);
if dimen==3
    s_11=zeros(params.nfreq+1,nobs);
    s_22=zeros(params.nfreq+1,nobs);
    s_33=zeros(params.nfreq+1,nobs);
    for p=1:params.nloop
        if(p>params.nwarmup)
            xi_curr=fit(fit(1).nexp_curr(p)).xi(:,p);
            spec_hat_curr_11=squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,:,p));
            spec_hat_curr_22=squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,:,p));
            spec_hat_curr_33=squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,3,:,:,p));
            for j=1:fit(1).nexp_curr(p)
                if(j==1)
                    s_11(:,1:xi_curr)=s_11(:,1:xi_curr(j))+repmat(spec_hat_curr_11(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                    s_22(:,1:xi_curr)=s_22(:,1:xi_curr(j))+repmat(spec_hat_curr_22(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                    s_33(:,1:xi_curr)=s_33(:,1:xi_curr(j))+repmat(spec_hat_curr_33(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                else
                    s_11(:,xi_curr(j-1)+1:xi_curr(j))=s_11(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_11(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                    s_22(:,xi_curr(j-1)+1:xi_curr(j))=s_22(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_22(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                    s_33(:,xi_curr(j-1)+1:xi_curr(j))=s_33(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_33(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                end
            end
        end
    end
    
    coh_21=zeros(params.nfreq+1,nobs);
    coh_31=zeros(params.nfreq+1,nobs);
    coh_32=zeros(params.nfreq+1,nobs);
    for p=1:params.nloop
        if(p>params.nwarmup)
            xi_curr=fit(fit(1).nexp_curr(p)).xi(:,p);
            spec_hat_curr_21=abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,1,:,:,p))).^2./...
            (squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,:,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,:,p)));
            spec_hat_curr_31=abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,1,:,:,p))).^2./...
                (squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,:,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,3,:,:,p)));
            spec_hat_curr_32=abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,2,:,:,p))).^2./...
                (squeeze(spect_matrices{fit(1).nexp_curr(p)}(3,3,:,:,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,:,p)));
            for j=1:fit(1).nexp_curr(p)
                if(j==1)
                    coh_21(:,1:xi_curr(j))=coh_21(:,1:xi_curr(j))+repmat(spec_hat_curr_21(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                    coh_31(:,1:xi_curr(j))=coh_31(:,1:xi_curr(j))+repmat(spec_hat_curr_31(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                    coh_32(:,1:xi_curr(j))=coh_32(:,1:xi_curr(j))+repmat(spec_hat_curr_32(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                else
                    coh_21(:,xi_curr(j-1)+1:xi_curr(j))=coh_21(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_21(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                    coh_31(:,xi_curr(j-1)+1:xi_curr(j))=coh_31(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_31(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                    coh_32(:,xi_curr(j-1)+1:xi_curr(j))=coh_32(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_32(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                end
            end
        end
    end   
    spectra{1}=s_11;spectra{2}=s_22;spectra{3}=s_33;
    coh{1}=coh_21;coh{2}=coh_31; coh{3}=coh_32;
elseif dimen==2
    s_11=zeros(params.nfreq+1,nobs);
    s_22=zeros(params.nfreq+1,nobs);
    for p=1:params.nloop
        if(p>params.nwarmup)
            xi_curr=fit(fit(1).nexp_curr(p)).xi(:,p);
            spec_hat_curr_11=squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,:,p));
            spec_hat_curr_22=squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,:,p));
            for j=1:fit(1).nexp_curr(p)
                if(j==1)
                    s_11(:,1:xi_curr)=s_11(:,1:xi_curr(j))+repmat(spec_hat_curr_11(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                    s_22(:,1:xi_curr)=s_22(:,1:xi_curr(j))+repmat(spec_hat_curr_22(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                else
                    s_11(:,xi_curr(j-1)+1:xi_curr(j))=s_11(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_11(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                    s_22(:,xi_curr(j-1)+1:xi_curr(j))=s_22(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_22(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                end
            end
        end
    end
    coh_21=zeros(params.nfreq+1,nobs);
    for p=1:params.nloop
        if(p>params.nwarmup)
            xi_curr=fit(fit(1).nexp_curr(p)).xi(:,p);
            spec_hat_curr_21=abs(squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,1,:,:,p))).^2./...
            (squeeze(spect_matrices{fit(1).nexp_curr(p)}(1,1,:,:,p)).*squeeze(spect_matrices{fit(1).nexp_curr(p)}(2,2,:,:,p)));
            for j=1:fit(1).nexp_curr(p)
                if(j==1)
                    coh_21(:,1:xi_curr(j))=coh_21(:,1:xi_curr(j))+repmat(spec_hat_curr_21(:,j),1,xi_curr(j))/(params.nloop-params.nwarmup);
                else
                    coh_21(:,xi_curr(j-1)+1:xi_curr(j))=coh_21(:,xi_curr(j-1)+1:xi_curr(j))+...
                        repmat(spec_hat_curr_21(:,j),1,xi_curr(j)-xi_curr(j-1))/(params.nloop-params.nwarmup);
                end
            end
        end
    end
   
    spectra{1}=s_11;spectra{2}=s_22;
    coh{1}=coh_21;  
end
