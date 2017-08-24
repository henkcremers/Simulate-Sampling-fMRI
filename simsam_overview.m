% =========================================================================
%% Overview Simulate Sampling and Statistical Power fMRI analyses 
% =========================================================================
% Dependencies: 
% https://nl.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
% http://www.fil.ion.ucl.ac.uk/spm/software/spm12/

% close all; clc; clear all

%% (1) load the set-up information.(Cluster location, size etc) 
% This is info for the weak diffuse scenario, load setupinfoSL.mat for the SL scenario. 
load setupinfoWD.mat
% load setupinfoSL.mat
 
%% (2) Set-up the data
setup  = simsam_setup_data(setupinfo);

%% (3) Generate the full sample data
simsam = simsam_generate_data(setup);

%% (4) Run the sampling stats.
[samplestats,sampledata] = simsam_samplestats(simsam);
%save samplestats.mat samplestats  
  
%% (5) Create and save some figures
simsam_saveplots(simsam,sampledata,setup)
