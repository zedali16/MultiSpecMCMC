function [param] = OptsMultiSpect(varargin)
% This function sets the optional input arguments for the function MCBSpec().
%
%  (I) ONLY DEFAULT PARAEMTERS
%       If only defult values for the parameters are desired, then either:
%
%           a) the SECOND argument in MultiSpec() can be left missing, or
%
%           b) params=setOptions() can be defined and used as the second
%           argument of MultiSpec().
%
% (II) USING NONDEFAULT PARAEMTERS
%       If some options other than the default are desired:
%
%           1) Set all default parameters using (Ia) above.
%
%           2) Change desired parameters.  
%
%
% PARAMETERS
%
%       nloop           -   The number of iterations run.
%                               Default: 6000.
%       nwarmup         -   The length of the burn-in.
%                               Default: 200.
%       nexp_max        -   The maximum number of segments allowed, adjust when
%                           time series get longer.
%                               Default: 10    
%       tmin            -   The minimum number of observation in each segment.
%                               Default:  100
%       prob_mml        -   The probability of small jump, 1-prob of big jump.
%                               Default: 0.8
%       nbasis          -   The number of linear smoothing spline basis functions.
%                               Default: 10
%       tau_up_limit    -   The normal variance prior for smoothing parameters.
%                               Default:  10^4.
%       sigmasqalpha    -   The smoothing parameter for alpha.
%                               Default: 10^4
%       init            -   Initial number of partitions
%                               Default: 3
%       nfreq           -   The number of frequencies for spectrm
%                               Default: 50
%       verb            -   An indicator if the iteration number should be printed
%                           at the completion of each interation.   1 is yes, 0 is
%                           no.
%       covdiag         -   An indicator if the diagnostic should be
%       ee              -   Step size in Hamiltonian Monte Carlo 
%                           Default is 0.1




param = struct('nloop',10000, 'nwarmup',2000, ...
               'nexp_max',10, 'tmin',60, 'prob_mml',0.8 , 'nbasis',10, ...
               'tau_up_limit',10^5 , 'sigmasqalpha',10^5, 'init',3,...
               'nfreq',50, 'verb',1, 'convdiag',0, 'ee',0.1);
               
param.bands = {};
               


 end
               