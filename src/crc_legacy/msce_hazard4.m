function val = msce_hazard4(model, t)
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

fun = @(s) ((q-p) ./  (q*exp(-p*s) - p*exp(-q*s))).^(r/a) - 1;

val = zeros(size(t));
val(1) = integral(fun, 0, t(1));
for i = 2:length(t)
    val(i) = val(i-1) + integral(fun, t(i-1), t(i));
end

val = - model.mu*expm1(model.mu*val);
end

