% Common Parameters
T = 100;
tt = linspace(0, T, 1001);
model = crc_mkmodel();
TOL = 1e-15;

%% Test Mk Survival 1 (normal, log)
S = crc_survival1(model, tt);
s = crc_mksurvival1(model, T);
logs = crc_mklogsurvival1(model, T);

assert(norm(S-s(tt), 'inf') < TOL);
assert(norm(S-exp(logs(tt)), 'inf') < TOL);


%% Test Mk Survival 2 (normal, minus 1, log)
S = crc_survival2(model, tt);
s = crc_mksurvival2(model, T);
sm1 = crc_mksurvival2m1(model, T);
logs = crc_mklogsurvival2(model, T);

assert(norm(S-s(tt), 'inf') < TOL);
assert(norm(S-(sm1(tt)+1), 'inf') < TOL);
assert(norm(S-exp(logs(tt)), 'inf') < TOL);


%% Test Mk Survival 3 (normal, minus 1, log)
S = crc_survival3(model, tt);
s = crc_mksurvival3(model, T);
sm1 = crc_mksurvival3m1(model, T);
logs = crc_mklogsurvival3(model, T);

assert(norm(S-s(tt), 'inf') < TOL);
assert(norm(S-(sm1(tt)+1), 'inf') < TOL);
assert(norm(S-exp(logs(tt)), 'inf') < TOL);


%% Test Mk Survival 4 (normal, log)
s = crc_mksurvival4(model, T);
logs = crc_mklogsurvival4(model, T);

assert(norm(s(tt)-exp(logs(tt)), 'inf') < TOL);