function times = poissproc_hom(lambda, T)
%POISSPROC_HOM Homogeneous Poisson process arrival times
%
% Generate the arrival times of a homogeneous poisson process
%
%   times = POISSPROC_HOM(lambda, T)
%
%   INPUT:
%       lambda: rate of process > 0, scalar
%       T: end time of process > 0, scalar
%   OUTPUT:
%       times: arrival times. >= 0, < T, vector
%
% See also poissproc
%
% Copyright 2018 Lukas Koestler (TUM) 

N = poissrnd(lambda*T);
times = sort(rand(1, N)*T);
end


