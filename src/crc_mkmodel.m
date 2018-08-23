function model = crc_mkmodel(alpha, X, p, q, mu, roa)
%CRC_MKMODEL Make a colorectal cancer model
%
% Model parameters for the crc model. See eqn. (5.6) in [1]
%
%   model = CRC_MKMODEL(alpha, X, p, q, mu, roa)
%
%   INPUT:
%       alpha:  as in [1]
%       X:      as in [1]
%       p:      as in [1]
%       q:      as in [1]
%       mu:     as in [1] (mu = mu_0 = mu_1)
%       roa:    as in [1] (rho over alpha; rho/alpha)
%   OUTPUT:
%       model: struct
%
%   DESCRIPTION:
%       If no parameters are given, i.e. model = CRC_MKMODEL(), then the
%       parameters from [1] are used
%
% See also crc_simulation
%
% Copyright 2018 Lukas Koestler (TUM)
%
% [1] Jihyoun Jeon et al. “Evaluation of screening strategies for
%     pre-malignant lesions using a biomathematical approach”.
%     In: Mathematical biosciences 213.1 (2008), pp. 56–70.

if nargin == 0
    alpha = 9;
    X = 1e8;
    p = -1.519930e-1;
    q = 3.893446e-6;
    mu = 1.364459e-6;
    roa = 6.886327;
end

% generate struct
model = struct(...
    'X', X,...
    'alpha', alpha,...
    'mu', mu,...
    'p', p, 'q', q,...
    'roa', roa);
end

