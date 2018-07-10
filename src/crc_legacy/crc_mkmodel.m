function model = crc_mkmodel()
%CRC_MKMODEL Make a colorectal cancer model

% Using parameters from paper
alpha = 9;
X = 1e8;
p = -1.519930e-1;
q = 3.893446e-6;
mu = 1.364459e-6;
rho = alpha*6.886327;

% generate struct
model = struct(...
    'X', X,...
    'alpha', alpha,...
    'mu', mu,...
    'p', p, 'q', q,...
    'rho', rho);
end

