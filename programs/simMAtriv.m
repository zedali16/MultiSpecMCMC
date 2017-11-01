function  TS0 = simMAtriv(T0, cut, ma11, ma12, ma21, ma22, Sigma, burn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a trivariate VMA time series
%
%   Input:
%       1) T - an integer >2
%       2) Phi1 - a P dim matrix for lag1
%       3) Phi2 - a P dim matrix for lag2
%       4) Sigma - PxP covariance matrix for error
%       5) burn - the size of the burn-in
%   Output:
%       1) 1xT matrix of the time series

T=T0+burn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save some info
error = mvnrnd([0;0;0],Sigma,T+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate starting terms

TS(:,1) = error(3,:) + (ma11*error(2,:)')' + (ma12*error(1,:)')';
TS(:,2) = error(4,:) + (ma11*error(3,:)')' + (ma12*error(2,:)')';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get the rest
for t = 3:T;
    if t<=cut+burn
        TS(:,t) = error(t+2,:) + (ma11*error(t+1,:)')' + (ma12*error(t,:)')';
    else
        TS(:,t) = error(t+2,:) + (ma21*error(t+1,:)')' + (ma22*error(t,:)')';
    end
end
TS0 = TS(:,(burn+1):end);
