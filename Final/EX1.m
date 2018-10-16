%% Ex1

clear all;
close all;
clc;
%% Load Dataset
load('ecoli_core_model.mat');

%% 1
% Aerobic Conditions
% Choose the solver
changeCobraSolver('mosek');

% Define the biomass reaction as the objective function.
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%% 2
% Set the max import rate of O2
model = changeRxnBounds(model,'EX_o2(e)',-20,'l');

%% 3&4
% Find all carbon sources
CarbonSources_vec = {'EX_glc(e)', 'EX_lac-D(e)', 'EX_fru(e)', 'EX_etoh(e)'};

% Set exchange rates of carbon sources to zero
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_lac-D(e)',0,'l');
model = changeRxnBounds(model,'EX_fru(e)',0,'l');
model = changeRxnBounds(model,'EX_etoh(e)',0,'l');


CarbonSources_Growth = zeros(1,4);

for i=1:4 % For all 4 carbon sources:
    % Set the import boundry to 10
    model = changeRxnBounds(model,CarbonSources_vec(i),-10,'l');
    
    % Get max biomass production by optimizing the model
    FBAsolution =  optimizeCbModel(model,'max');
    
    % Find max growth rate
    CarbonSources_Growth(i)= FBAsolution.f;
    
    % Reset the boundary of the source to 0
    model = changeRxnBounds(model,CarbonSources_vec(i),0,'l');
   
end

%% 5
% Anaerobic Conditions
% Set the import rate of O2 to 0
model = changeRxnBounds(model,'EX_o2(e)',0,'l');


CarbonSources_Growth_A = zeros(1,4);

% Set exchange rates of carbon sources to zero
model = changeRxnBounds(model,'EX_glc(e)',0,'l');
model = changeRxnBounds(model,'EX_lac-D(e)',0,'l');
model = changeRxnBounds(model,'EX_fru(e)',0,'l');
model = changeRxnBounds(model,'EX_etoh(e)',0,'l');

for i=1:4 % For all 4 carbon sources:
    % Set the import boundry to 10
    model = changeRxnBounds(model,CarbonSources_vec(i),-10,'l');
    
    % Get max biomass production by optimizing the model
    FBAsolution =  optimizeCbModel(model,'max');
    
    % Find max growth rate
    CarbonSources_Growth_A(i)= FBAsolution.f;
    
    % Reset the boundary of the source to 0
    model = changeRxnBounds(model,CarbonSources_vec(i),0,'l');
   
end

%% 6
% This step is done using the Escher map tool on github. The results are
% reported in the group report

%% 7
% The filled table is presented in the report