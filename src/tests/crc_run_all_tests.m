curr_filename = mfilename();
curr_fullpath = mfilename('fullpath');
curr_folderpath = curr_fullpath(1:end-length(curr_filename));
fprintf('Running all tests in: %s\n\n',curr_folderpath);

suite = matlab.unittest.TestSuite.fromFolder(curr_folderpath);
result = run(suite)