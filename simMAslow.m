function  TS0 = simMAslow(T0, cont, Sigma, burn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/71/2011 for the Peanalized Multivariate Whittle Liklihood paper
% Generate a VAR TS
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
a1=cont*(1-1.781*sin(pi*(1:T0)/(2*T0)));
a2=cont*(1-1.781*cos(0.8*pi*(1:T0)/(2*T0)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save some info
P = size(Sigma,1);
error=mvnrnd([0;0],Sigma,T+2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate starting terms
a=[a1(1) -1; -1, a2(1)];
TS(:,1) = error(3,:)' + a*error(2,:)' + diag([0.5 -1.2])*error(1,:)';
a=[a1(2) -1; -1, a2(2)];
TS(:,2) = error(4,:)' + a*error(3,:)' + diag([0.5 -1.2])*error(2,:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get the rest
for t = 3:T
    a=[a1(t) -1; -1 a2(t)];
    TS(:,t) = error(t+2,:)' + a*error(t+1,:)' + diag([0.5 -1.2])*error(t,:)';
end

TS0 = TS(:,(burn+1):end);