function fun = crc_mklogsurvival1(model, T)
%CRC_MKLOGSURVIVAL1 log(S_1(t)) function ([1] eqn. (C.3))
%
%   INPUTS:
%       model: parameters from tsce_mkmodel, struct
%   OUTPUT:
%       fun: survival function as chebfun for t in [0, 100]
%
% See also crc_mkmodel
%
% [1]   Jihyoun Jeon et al. “Evaluation of screening strategies for
%       pre-malignant lesions using a biomathematical approach”.
%       In: Mathematical biosciences 213.1 (2008), pp. 56–70.

if nargin == 0
    model = crc_mkmodel();
end

if nargin <= 1
    T = 50;
end
    


loc_fun = @(t) (model.p*model.q)/(model.alpha) * (exp(-model.p*t) - exp(-model.q*t)) ./ ...
    (model.q*exp(-model.p*t) - model.p*exp(-model.q*t));
fun = log1p(chebfun(loc_fun, [0, T]));
end

