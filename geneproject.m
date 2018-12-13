clear all;
close all;
clc;

%% Load Dataset and Set up model
load('photosynthesis.mat');

% Choose the solver
changeCobraSolver('mosek');

% Define the biomass reaction as the objective function.
model = changeObjective(model,'Ec_biomass_SynHetero');

%% Wild type
% The network viability criterion was a minimum growth rate of 10% of the wild-type growth. 

wildSol = optimizeCbModel(model,'max');
wildGrowth = wildSol.f;

%% autotrophic conditions
%  For autotrophic growth,the carbon dioxide uptake was set to 3.7 mmol/gDW/h and the photon uptake to 100 mmol/gDW/h.

% for i=1:50
%     
%     model1 = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
%     model1 = changeRxnBounds (model1, 'EX_glc(e)', 0,'l');
%     model1 = changeRxnBounds (model1, 'EX_hco3(e)',-3.7,'l');
%
%     flag = 0;
%     num_essenGenes1=0;
%     now_genes = model.genes;
%     
%     % simulation
%     while flag<50;
%                 
%         rng(0,'twister');   
%         gindex = randi([1,length(now_genes)]);
%         deleteGene = now_genes(gindex);
%         
%         %delete sol: [model, hasEffect, constrRxnNames, deletedGenes] 
%         deleteSol = deleteModelGenes(model1, deleteGene);
%         isol = optimizeCbModel(deleteSol.model);
%         
%         if isol.f >= 0.9 * wildSol.f
%             flag = 0;
%         end
%         
%         if isol.f < 0.9 * wildSol.f
%             model1 = deleteSol.model;
%             now_genes(gindex)=[];
%             flag ++;
%         end
%     end
%        
%     num_essential(i)=length(now_genes);    
% end
    
     model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
     model = changeRxnBounds (model, 'EX_glc(e)', 0,'l');
     model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');
     
     Sol1 = optimizeCbModel(model,'max');
     

     [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
     
     num1=0;
     EssentialGene1{1}= 'autotrophic';
     

    for j=1:length(grRatio)
      if grRatio(j)<0.1
          num1=num1+1;
          EssentialGene1{num1+1}= model.genes{j};
      end
    end
   

%% mixotrophic condition
% For mixotrophic growth, glucose was added to the system by allowing an uptake corresponding to 0.85 mmol/gDW/h. 

model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -0.85,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');

Sol2 = optimizeCbModel(model,'max');
     

[grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
     
num2=0;
EssentialGene2{1}= 'mixotrophic';

    for j=1:length(grRatio)
      if grRatio(j)<0.1
          num2=num2+1;
          EssentialGene2{num2+1}= model.genes{j};
      end
    end


%% heterotrophic condition
%  For heterotrophic growth, light uptake was set to zero and glucose uptake to 0.85 mmol/gDW/h. 

model = changeRxnBounds (model, 'EX_photon(e)', 0,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -0.85,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');

Sol3 = optimizeCbModel(model,'max');
     

[grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
     
num3=0;
EssentialGene3{1}= 'heterotrophic';

    for j=1:length(grRatio)
      if grRatio(j)<0.1
          num3=num3+1;
          EssentialGene3{num3+1}= model.genes{j};
      end
    end



%% Simulation 
% A single simulation starts with a full network and selects a gene at random.
% If removal of this particular gene results in a viable network, the gene is permanently removed. 
% The process is repeated for all the genes in the network. 
% The remaining genes are essential for network survival. 
% On the basis of 1000 simulation runs, sets of always-present genes (essential), sometimes-present genes (synthetic lethal) and never-present (nonessential) were identified.

