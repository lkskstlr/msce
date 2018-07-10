function fun = crc_mksurvival3(model, T)
%CRC_MKSURVIVAL3 S_3(t) function ([1] eqn. (C.1) with k=3)
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

logS3 = crc_mklogsurvival3(model, T);
fun = chebfun(@(t) exp(logS3(t)), [0, T]);
end

