function samples = helper_sample_h(h, sigma, epsilon, num_samples)
%HELPER_SAMPLE_H Sample a very specific distribution
%
%   samples = HELPER_SAMPLE_H(h, sigma, epsilon, num_samples)
%
%   INPUT:
%       h: cdf that should be sampled, chebfun, h(0) = 0, h(sigma) = 1
%       sigma: 0 <= samples <= sigma, upper bound of interval
%       epsilon: for 0<=x<=epsilon an approximation is used
%               choose as small as possible as long the run-time doesn't
%               explode. 0.001*sigma might be a good value
%       num_samples: number of samples from h
%   OUTPUT:
%       samples: the samples
%
%   DESCRIPTION:
%           h is the mean value function of the two-stage process
%           it is assumed that h = (eta * lambda) (the * means convolution)
%           where the function lambda is essentially constant (of Order 1)
%           and the function eta is the integral over an essentially
%           constant function, i.e. essentially linear. Thus h is
%           essentially a parabola. Because h resembles the cdf, the
%           underlying pdf is essentially a linear function. These
%           assumptions are used HEAVILY and using this function for
%           another purpose will fail!
%
%
% Copyright 2018 Lukas Koestler (TUM)


%% build approximation for 0 <= x <= epsilon
% h(0) = 0, h'(0) = 0
% fit h(epsilon) and h'(epsilon)
% happrox(x) = a x^4 + b x^2 is easy to invert with z = x^2
% happrox(x) = a x^3 + b x^2 would be harder to invert

dh = diff(h);
A = [epsilon^4, epsilon^2;4*epsilon^3, 2*epsilon];
hepsilon = h(epsilon);
rhs = [hepsilon;dh(epsilon)];
coeff = A \ rhs;

a = coeff(1);
b = coeff(2);
if (a <= 0) || (b <= 0)
    error("Assumptions violated");
end

hsqrtinvapprox = @(y) sqrt((-b + sqrt(b^2 + 4*a*y.^2)) / (2*a));

%% invert sqrt(h) on epsilon < x <= sigma
% This is necessary because inverting sqrt(h) around x = 0 is hard
% because h'(0) = 0. Although the sqrt(.) aleviates the effect
% (which is why I do it this way) on has to use an approximation
% for small x: x <= epsilon
hsqrtinv = inv(sqrt(restrict(h, [epsilon, sigma])));

%% Sample
samples = zeros(1, num_samples);

% sample linear distribution for u, i.e p(u) = 2*u
% this helps because h' approx (t) and thus the distrbutions fit well.
% Therefore I need inv(sqrt(h)) and not inv(h) as for the normal
% transform
u = sqrt(rand(1, num_samples));

ind_approx = u <= sqrt(hepsilon);
samples(ind_approx) = hsqrtinvapprox(u(ind_approx));
samples(~ind_approx) = hsqrtinv(u(~ind_approx));
end

