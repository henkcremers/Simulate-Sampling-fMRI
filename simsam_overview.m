% =========================================================================
%% Overview Simulate Sampling and Statistical Power fMRI analyses 
% =========================================================================
% Dependencies: 
% https://nl.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
% http://www.fil.ion.ucl.ac.uk/spm/software/spm12/

% close all; clc; clear all

%% (1) load the set-up information.(Cluster location, size etc) 
% this is an example of the weak diffuse scenario. Alternatively you can
% generate a new set of brain-behavior correlations and main effects 
% (see simsam_setup_scenario.m for examples). 
load setupinfoWD.mat
 
%% (2) Set-up the data
setup  = simsam_setup_data(setupinfo);

%% (3) Generate the full sample data
simsam = simsam_generate_data(setup);

%% (4) Run the sampling stats.
[samplestats,sampledata] = simsam_samplestats(simsam);
%save samplestats.mat samplestats  
  
%% (5) Create and save some figures
simsam_saveplots(simsam,sampledata,setup)
