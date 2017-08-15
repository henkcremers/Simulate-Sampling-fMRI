%%  Set up clusters for different Brain-Behavior correlation scenarios

% orientation (hor slice/#35): left-right: 10-70. down-up: 10-90. 

%% Brain-Behavior Correlation - r
% Weak Diffuse Scenario
% cl1 = [16 30 10 1];
% cl2 = [15 70 5 1; 15 85 5 1; 18 65 8 1; 20 90 4 1]; 
% cl3 = [35 20 10 -1; 45 23 10 -1; 60 30 10 -1; 50 45 5 -1; 45 40 8 -1];
% cl4 = [58 80 12 -1; 62 66 10 -1; 60 60 6 -1]; 
% cl5 = [30 88 15 1; 40 75 15 1; 45 90 15 1 ];
% cl6 = [10 30 12 1; 30 40 12 1; 8 28 22 1];
% cl7 = [10 50 14 -1; 10 60 15 -1];
% cl8 = [70 45 12 1];
% clusters = [cl1;cl2;cl3;cl4;cl5;cl6;cl7;cl8];
% range = [0.05 0.15];
% smooth = 4;

% setupinfo.BBclusters = clusters;
% setupinfo.BBsmooth = 4;
% setupinfo.BBrange = [0.05 0.15];
% setupinfo.BBscenario = 'WeakDiffuse';
% setupinfo.orientation = 'hor';
% setupinfo.slice = 35;

% Strong Localised
% ----------------
range = [0.7 0.9];
cl2 = [16 71 5 1];
cl5 = [55 28 5 -1];
cl6 = [38 90 5 1];
cl7 = [16 30 5 -1];
clusters = [cl2;cl5;cl6;cl7];
scenario = 'StrongLocalised';
smooth = 2;

setupinfo.BBclusters = clusters;
setupinfo.BBsmooth = smooth;
setupinfo.BBrange = range;
setupinfo.BBscenario = scenario;
setupinfo.orientation = 'hor';
setupinfo.slice = 35;

%% Main-effects - Cohen's D
% ---------------------------------------------------------------
cl1 = [18 75 7 1; 60 65 7 1]; 
cl2 = [38 35 7 -1; 42 30 6 -1];
clusters = [cl1;cl2];

setupinfo.MEclusters = clusters;
setupinfo.MEsmooth = 3;
setupinfo.MErange = [0.2 0.5];

% save the info
%mkdir(setupinfo.BBscenario)
%cd(setupinfo.BBscenario);
save setupinfoSL.mat setupinfo
