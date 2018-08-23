function val = crc_hazard3(model, t)
%CRC_HAZARD3 h_3(t) function ([1] eqn. (C.7))
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


num = q-p;
den = q*exp(-p*t) - p*exp(-q*t);

val = model.mu*(1 - (num./den).^(model.roa));
end

