%% Ex2

clear all;
close all;
clc;

%% Load Dataset
load('ecoli_core_model.mat');

%% Pre-set
%Choose the solver
changeCobraSolver('mosek')

% Define the biomass reaction as the objective function.
model= changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%Define the subtrates to be analyzed
Sub={'EX_glc(e)','EX_lac-D(e)','EX_fru(e)','EX_etoh(e)','EX_gln-L(e)','EX_glu-L(e)','EX_pyr(e)','EX_ac(e)'}; 

%Create the set for data storage
for n=1:8
    EssentialGenes{1,2*n-1}=Sub{n};
    EssentialGenes{1,2*n}=Sub{n};
    EssentialGenes{2,2*n-1}='Aerobic';
    EssentialGenes{2,2*n}='Anaerobic';
end

%% 1&2

%Analyze each subtrate
for n=1:8
    
    % Set the max rate for previous subtrate
    if n~=1
        model=changeRxnBounds(model,Sub{n-1},0,'l');  
    end
    
    % Set the max rate of current subtrate
    model= changeRxnBounds(model,Sub{n},-10,'l');

    % Set the gene knockout function (Aerobic)
    model = changeRxnBounds(model,'EX_o2(e)',-20,'l');
    [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
    
    %Analyze and store the essential genes and total number for this subtrate
    num_aerobic=0;
    for i=1:length(grRatio)
       if grRatio(i)<0.1
           EssentialGenes{num_aerobic+4,2*n-1}= model.genes{i};
           num_aerobic=num_aerobic+1; 
        end
    end
    EssentialGenes{3,2*n-1}=num_aerobic;
    
    % Set the gene knockout function (Anaerobic)
    model = changeRxnBounds(model,'EX_o2(e)',0,'l');
    [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
     
    % Analyze and store the essential genes and total number for this subtrate
    num_anaerobic=0;
    for i=1:length(grRatio)
      if grRatio(i)<0.1
          EssentialGenes{num_anaerobic+4,2*n}= model.genes{i};
          num_anaerobic=num_anaerobic+1;
       end
    end
    EssentialGenes{3,2*n}=num_anaerobic;
end   

%% 3 

%Create a table for the results and present
Table=table(EssentialGenes);
open Table