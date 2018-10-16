%% Ex3

clear all;
close all;
clc;

%% Load Dataset
load('ecoli_core_model.mat');

%% Pre-set
% Choose the solver
changeCobraSolver('mosek')

% Define the biomass reaction as the objective function.
model= changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%% 1 Aerobic
%Set the rate for O2 and glucose (Aerobic)
model = changeRxnBounds(model,'EX_o2(e)',-1000,'l');
model= changeRxnBounds(model,'EX_glc(e)',-10,'l');
 

%optimal growth is set to be 80% of the maximal growth yield
FBAsolution = optimizeCbModel(model,'max');
[minFlux, maxFlux] = fluxVariability(model, 80);

% Store the results
gluc_Aer_FluxVar=table(model.rxns,minFlux,maxFlux);

%% 2 Anaerobic
%Set the rate for O2 (Anaerobic)
model = changeRxnBounds(model,'EX_o2(e)',0,'l');

%optimal growth is set to be 80% of the maximal growth yield
FBAsolution = optimizeCbModel(model,'max');
[minFlux, maxFlux] = fluxVariability(model, 80);

% Store the results
gluc_anAer_FluxVar=table(model.rxns,minFlux,maxFlux);

%% 3
open gluc_anAer_FluxVar



