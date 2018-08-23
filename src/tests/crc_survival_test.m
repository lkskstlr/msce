% Common Parameters
model = crc_mkmodel();
load("crc_survival_test.mat");

%% Test Survival1: t = 0
assert(isequal(crc_survival1(model, 0), 1.0));

%% Test Survival1: t -> infty
assert(isequal(crc_survival1(model, 1000), (1+model.p/model.alpha)));

%% Test Survival1: t in [0, 100]
assert(isequal(crc_survival1(model, tt), S1));



%% Test Survival2: t = 0
assert(isequal(crc_survival2(model, 0), 1.0));

%% Test Survival2: t -> infty
assert(isequal(crc_survival2(model, 10000), 0));

%% Test Survival2: t in [0, 100]
assert(isequal(crc_survival2(model, tt), S2));



%% Test Survival3: t = 0
assert(isequal(crc_survival3(model, 0), 1.0));

%% Test Survival3: t -> infty
assert(isequal(crc_survival3(model, 1000000000), 0));

%% Test Survival3: t in [0, 100]
assert(isequal(crc_survival3(model, tt), S3));