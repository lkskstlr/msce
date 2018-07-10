function val = msce_survival1(model, t)
%MSCE_SURVIVAL1 S_1(t) function ([1] eqn. (C.3))
%
%   INPUTS:
%       model: parameters from tsce_mkmodel, struct
%       t: time, can be any array
%   OUTPUT:
%       val: value of survival function, same size as t
%
% See also tsce_mkmodel
%
% [1]   Jihyoun Jeon et al. “Evaluation of screening strategies for
%       pre-malignant lesions using a biomathematical approach”.
%       In: Mathematical biosciences 213.1 (2008), pp. 56–70.

num = exp(-model.p*t) - exp(-model.q*t);
den = model.q*exp(-model.p*t) - model.p*exp(-model.q*t);
val = 1 + (model.p*model.q*num)./(model.alpha*den);
end

