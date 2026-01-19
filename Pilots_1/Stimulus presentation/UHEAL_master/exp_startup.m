% startup
rootdir = fileparts(which('exp_startup.m'));
datadir = ([rootdir,'/_data']);


try %#ok
    rng(1); 
end  


fprintf('\n project directory now added to the current path \n')
addpath(genpath('_exp'))
addpath(genpath('_demo'))

fprintf('\n directory addded to the path')



