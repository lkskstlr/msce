function fun = crc_mklogsurvival4(model, T)
%CRC_MKLOGSURVIVAL4 log(S_4(t)) function ([1] eqn. (C.1)) as chebfun
%
%   INPUTS:
%       model: parameters from crc_mkmodel, struct
%   OUTPUT:
%       fun: log survival function as chebfun for t in [0, T]
%
% See also crc_mkmodel, crc_mklogsurvival1, crc_mksurvival1,
%                       crc_mklogsurvival2, crc_mksurvival2
%                       crc_mklogsurvival3, crc_mksurvival3,
%                       crc_mklogsurvival3, crc_mksurvival4
%
% [1] Jihyoun Jeon et al. “Evaluation of screening strategies for
%     pre-malignant lesions using a biomathematical approach”.
%     In: Mathematical biosciences 213.1 (2008), pp. 56–70.

if nargin == 0
    model = crc_mkmodel();
end

if nargin <= 1
    T = 50;
end

logS3 = crc_mklogsurvival3(model, T);
S3m1 = expm1(logS3);
S3m1h = cumsum(S3m1);
fun = model.mu*S3m1h;
end

