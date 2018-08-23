function val = crc_hazard4(model, t)
%CRC_HAZARD4 h_4(t) function ([1] eqn. (C.8))
%
%   INPUTS:
%       model: parameters from crc_mkmodel, struct
%       t: time >= 0, can be any array
%   OUTPUT:
%       val: value of hazard function, same size as t
%
% See also crc_mkmodel, crc_hazard1, crc_hazard2, crc_hazard3, crc_hazard4
%
% [1] Jihyoun Jeon et al. “Evaluation of screening strategies for
%     pre-malignant lesions using a biomathematical approach”.
%     In: Mathematical biosciences 213.1 (2008), pp. 56–70.

p = model.p;
q = model.q;

fun = @(s) ((q-p) ./  (q*exp(-p*s) - p*exp(-q*s))).^(model.roa) - 1;

val = zeros(size(t));
val(1) = integral(fun, 0, t(1));
for i = 2:length(t)
    val(i) = val(i-1) + integral(fun, t(i-1), t(i));
end

val = - model.mu*expm1(model.mu*val);
end

