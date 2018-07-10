function fun = crc_mklogsurvival2(model, T)
%CRC_MKLOGSURVIVAL2 log(S_2(t)) function ([1] eqn. (C.2))
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

f1 = chebfun(@(t) log(model.q*exp(-model.p*t) - model.p*exp(-model.q*t)),...
    [0, T]);
fun = (model.mu/model.alpha)*(log(model.q-model.p) - f1);

end

