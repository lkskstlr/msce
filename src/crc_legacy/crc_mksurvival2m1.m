function fun = crc_mksurvival2m1(model, T)
%CRC_MKSURVIVAL2 S_2(t)-1 function ([1] eqn. (C.2))
%
%   INPUTS:
%       model: parameters from crc_mkmodel, struct
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

fun = expm1(crc_mklogsurvival2(model, T));
end

