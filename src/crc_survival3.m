function val = crc_survival3(model, t)
%CRC_SURVIVAL3 S_3(t) function ([1] eqn. (C.1) with k=3)
%
%   INPUTS:
%       model: parameters from crc_mkmodel, struct
%       t: time, can be any array
%   OUTPUT:
%       val: value of survival function, same size as t
%
% See also crc_mkmodel, crc_survival1, crc_survival2, crc_survival3
%
% [1] Jihyoun Jeon et al. “Evaluation of screening strategies for
%     pre-malignant lesions using a biomathematical approach”.
%     In: Mathematical biosciences 213.1 (2008), pp. 56–70.

% Tolerances for integral
RelTol = 1e-6;
AbsTol = 1e-10;

[t_sorted, ind] = sort(t);
inv_ind(ind) = 1:length(t_sorted);

S2_int = zeros(size(t_sorted));

S2_fun = @(s) crc_survival2(model, s);
S2_int(1) = integral(S2_fun, 0, t_sorted(1));

for i = 2:length(t_sorted)
    S2_int(i) = S2_int(i-1) + integral(S2_fun, t_sorted(i-1), t_sorted(i),...
        'RelTol', RelTol, 'AbsTol', AbsTol);
end

val = exp(model.mu*( S2_int(inv_ind) - t ));
end

