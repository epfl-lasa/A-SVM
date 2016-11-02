% Script to setup the matlab path for calling SVMWidget functions


disp('Adding functions to matlab PATH');
rt = pwd;
addpath(genpath([rt '/matlab']), '-begin');
disp('Done.');