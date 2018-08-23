%% Test Basic Functionality
model1 = crc_mkmodel();
model2 = crc_mkmodel(9, 1e8, -1.519930e-1,...
    3.893446e-6, 1.364459e-6, 6.886327);
assert(isequal(model1, model2))