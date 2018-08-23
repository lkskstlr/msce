function val = crc_hazard1(model, t)
%CRC_HAZARD1 h_1(t) function ([1] eqn. (C.5))
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

a = model.alpha;
p = model.p;
q = model.q;

et = exp((q-p)*t);

num = p*q*(p-q)^2*et;
den = (q*(a+p)*et - p*(a+q)).*(q*et-p);

val = - num ./ den;
end

