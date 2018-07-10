function val = msce_hazard2(model, t)
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

a = model.alpha;
p = model.p;
q = model.q;
r = model.rho;

et = exp((q-p)*t);

num = r*p*q*(et-1);
den = a*(q*et - p);

val = - num ./ den;
end

