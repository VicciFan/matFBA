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

% wildSol = optimizeCbModel(model,'max');
% wildGrowth = wildSol.f;


%% autotrophic conditions   HCO3  
model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
model = changeRxnBounds (model, 'EX_glc(e)', 0,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');
model = changeRxnBounds (model, 'EX_co2(e)',0,'l');

num = 0;

for i= -6:0.5:4
    model = changeRxnBounds (model, 'EX_h(e)',i,'l');
    %model = changeRxnBounds (model, 'EX_h(e)',i,4);
    Sol = optimizeCbModel(model,'max');
    num=num+1;
    RateHCO3(num) = Sol.f;
end

    
%% autotrophic conditions   CO2
model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
model = changeRxnBounds (model, 'EX_glc(e)', 0,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',0,'l');
model = changeRxnBounds (model, 'EX_co2(e)',-3.7,'l');


num = 0;

for i= -6:0.5:4
    model = changeRxnBounds (model, 'EX_h(e)',i,'l');
    Sol = optimizeCbModel(model,'max');
    num=num+1;
    RateCO2(num) = Sol.f;
end


%% heterotrophic condition 


model = changeRxnBounds (model, 'EX_photon(e)', 0,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -0.85,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',0,'l');
model = changeRxnBounds (model, 'EX_co2(e)',0,'l');

num = 0;

for i= -6:0.5:4
    model = changeRxnBounds (model, 'EX_h(e)',i,'l');
    Sol = optimizeCbModel(model,'max');
    num=num+1;
    RateHetero(num) = Sol.f;
end
%% mixotrophic condition

model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -0.38,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');
model = changeRxnBounds (model, 'EX_co2(e)',0,'l');

num = 0;

for i= -6:0.5:4
    model = changeRxnBounds (model, 'EX_h(e)',i,'l');
    Sol = optimizeCbModel(model,'max');
    num=num+1;
    RateMix(num) = Sol.f;
end


%%
x= linspace( -6,4,21);
plot (x,RateHCO3/max(RateHCO3));
hold on 
plot (x,RateCO2/max(RateCO2));
hold on 
plot (x,RateHetero/max(RateHetero));
hold on
plot (x,RateMix/max(RateMix));
legend('HCO3','CO2','Hetero','Mix')







%%
% The exchange of protons between the cell and the medium was varied from -6 to 4 mmol.gDW-1.h-1
% the relative growth rate was computed under heterotrophic(blue line), 
% mixotrophic (orange line) 
% and autotrophic conditions, using CO2 (green solid line) and HCO3 (green dots line) as inorganic carbon, respectively
% %% heterotrophic condition (Nitrogen)


%% 
% model = changeRxnBounds (model, 'EX_photon(e)', 0,'l');
% model = changeRxnBounds (model, 'EX_glc(e)', -1,'l');
% model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');
% 
% 
% model = changeRxnBounds (model, 'EX_no3(e)', 0,'l');
% for i=1:30
% 
%     model = changeRxnBounds (model,'EX_'+ Carbon{i},-1,'l');
%     
%     SolHeteroNi{i} = optimizeCbModel(model,'max');
%     
%     if SolHeteroNi{i}.f >= 0.1 * SolHeteroWild.f
%         SolHeteroNiYN{i} = SolHeteroNi{i}.f/SolHeteroWild.f ;
%     else
%         SolHeteroNiYN{i} = 0;
%     end
%     
%     model = changeRxnBounds (model,'EX_'+ Carbon{i},0,'l');
% end


