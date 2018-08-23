% Common Parameters
model = crc_mkmodel();
load("crc_hazard_test.mat");


%% Test Hazard 1: t in [0, 100]
assert(isequal(crc_hazard1(model, tt), h1));

%% Test Hazard 2: t in [0, 100]
assert(isequal(crc_hazard2(model, tt), h2));

%% Test Hazard 3: t in [0, 100]
assert(isequal(crc_hazard3(model, tt), h3));

%% Test Hazard 4: t in [0, 100]
assert(isequal(crc_hazard4(model, tt), h4));