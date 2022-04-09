function clph2Install
% add clph2 into the MATLAB search path

here = pwd;
addpath(here);
addpath([here,filesep,'utils']);
savepath

fprintf('\n clph2 has been added to MATLAB search path.\n');
fprintf(' run clph2Test to see performance.\n');

end
