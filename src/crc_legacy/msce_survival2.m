function val = msce_survival2(model, t)
%MSCE_SURVIVAL2 S_2(t) function ([1] eqn. (C.2))
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


num = model.q-model.p;
den = model.q*exp(-model.p*t) - model.p*exp(-model.q*t);
val = (num ./ den).^(model.mu/model.alpha); 
end

