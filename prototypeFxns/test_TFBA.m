% Test TFBA algorithm 18/07/18
clc,clearvars,close all

%% Example 1
% Simple pathway with 1 balanced (internal) metabolite, 3 internal fluxes
% and 3 non-balanced (external) metabolites
clc,clearvars,close all
S = [-1,0,0;...
    1,-1,-1;...
    0,1,0;...
    0,0,1];
intMets = boolean([0,1,0,0])';    % Only X2 is only met balanced
vmin    = [6,-2,0]';              % [mmol/gdc/h]
vmax    = [10,6,15]';             % [mmol/gdc/h]
xmin    = [1e-4,1e-5,1e-4,1e-6]'; % [M]
xmax    = [1e-3,1e-3,1e-3,1e-4]'; % [M]
DGr_std = [-15,-5,-10]';          % [kJ/mol] (*)
% (*): Note that if the rxns that exchange mass with the surroundings are 
% regarded as transport reactions, this change should be reflected on how
% the DGr_std is calculated (e.g., membrane potential has to be considered,
% refer to Henry et al. 2007 Biophys J)

% Run TFBA
[DGr_rng,xrng,vrng] = TFBA(S,DGr_std,vmin,vmax,xmin,xmax,intMets);