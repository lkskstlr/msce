function [N2, N3, N4, N_polyp] = crc_simulation(model, sigma, num_samples)
%CRC_SIMULATION
%
%   INPUTS:
%       model:          parameters from crc_mkmodel, struct
%       sigma:          age at screen
%       num_samples:    number of samples
%   OUTPUT:
%       N2:             Number of APC+/- cells per individual
%       N3:             Number of APC-/- cells per sample
%       N4:             Number of Polyp cells per sample
%       N_polyp:        cell array, N_polyp{i} = sizes of polyps for
%                       individual i
%
% See also crc_mkmodel
%
% [1]   Jihyoun Jeon et al. “Evaluation of screening strategies for
%       pre-malignant lesions using a biomathematical approach”.
%       In: Mathematical biosciences 213.1 (2008), pp. 56–70.

%% Chebfun
S3 = crc_mksurvival3(model, sigma);
S2 = crc_mksurvival2(model, sigma);
lambda3 = model.mu*model.X*chebfun(@(tau) S3(sigma-tau), [0, sigma]);
lambda2 = model.mu*chebfun(@(tau) S2(sigma-tau), [0, sigma]);

%% Simulate APC +/- and APC -/- cells
[N2, N3, ind, times_apc_mm] = poissproc2(...
    cumsum(lambda3),...
    lambda2,...
    sigma,...
    num_samples);

%% Simulate Polyps
% Simulate Polyps
% (3.30)    R = rho*X/alpha (here X = 1)
%           P = 1 - alpha*zeta
%           SizePDF = negative binomial (R,P)
numer = exp(-model.p*(sigma-times_apc_mm))...
    - exp(-model.q*(sigma-times_apc_mm));
denum =...
    (model.q+model.alpha)*exp(-model.p*(sigma-times_apc_mm))-...
    (model.p+model.alpha)*exp(-model.q*(sigma-times_apc_mm));
zeta =numer ./ denum;

PolypSizes =  nbinrnd(model.rho/model.alpha, 1 - model.alpha*zeta);

%% Transform to correct datastructure
N4 = zeros(1, num_samples);
N_polyp = cell(1, num_samples);

for i = 1:num_samples
    N_polyp{i} = PolypSizes(ind(i):ind(i+1)-1);
    N4(i) = sum(PolypSizes(ind(i):ind(i+1)-1));
end

end

