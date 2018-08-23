function val = crc_survival2(model, t)
%CRC_SURVIVAL2 S_2(t) function ([1] eqn. (C.2))
%
%   INPUTS:
%       model: parameters from crc_mkmodel, struct
%       t: time>=0, can be any array
%   OUTPUT:
%       val: value of survival function, same size as t
%
% See also crc_mkmodel, crc_survival1, crc_survival2, crc_survival3
%
% [1] Jihyoun Jeon et al. “Evaluation of screening strategies for
%     pre-malignant lesions using a biomathematical approach”.
%     In: Mathematical biosciences 213.1 (2008), pp. 56–70.


num = model.q-model.p;
den = model.q*exp(-model.p*t) - model.p*exp(-model.q*t);
val = (num ./ den).^(model.mu/model.alpha); 
end

