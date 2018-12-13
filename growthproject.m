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

%% Set up source

Carbon = {'fru[c]','ac[c]','fum[c]','succ[c]','cit[c]','akg[c]','mal-L[c]','pyr[c]','arg-L[c]','gly[c]','his-L[c]','ser-L[c]'};
    
%% heterotrophic condition carbon


model = changeRxnBounds (model, 'EX_photon(e)', 0,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -1,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');

SolHeteroWild = optimizeCbModel(model,'max');


model = changeRxnBounds (model, 'EX_glc(e)', 0,'l');

for i=1:length(Carbon)
    met = Carbon{i};
    model = addSinkReactions(model, {met},0,'l');
    model = changeRxnBounds (model,['sink_',met],-1,'l');
    
    SolHeteroCar{i} = optimizeCbModel(model,'max');
    if SolHeteroCar{i}.f >= 0.1 * SolHeteroWild.f
        SolHeteroCarYN{i} = SolHeteroCar{i}.f/SolHeteroWild.f ;
    else
        SolHeteroCarYN{i} = 0;
    end
    
    model = changeRxnBounds (model,['sink_',met],0,'l');
end


%% mix carbon

model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -1,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');

SolMixWild = optimizeCbModel(model,'max');


model = changeRxnBounds (model, 'EX_glc(e)', 0,'l');

for i=1:length(Carbon)
    met = Carbon{i};

    model = changeRxnBounds (model,['sink_',met],-1,'l');
    
    SolMixCar{i} = optimizeCbModel(model,'max');
    if SolMixCar{i}.f >= 0.1 * SolHeteroWild.f
        SolMixCarYN{i} = SolMixCar{i}.f/SolMixWild.f ;
    else
        SolMixCarYN{i} = 0;
    end
    
    model = changeRxnBounds (model,['sink_',met],0,'l');
end

%% heterotrophic condition (Nitrogen)


model = changeRxnBounds (model, 'EX_photon(e)', 0,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -1,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');


model = changeRxnBounds (model, 'EX_no3(e)', 0,'l');

for i=1:length(Carbon)
    met = Carbon{i};
    model = changeRxnBounds (model,['sink_',met],-1,'l');
    
    SolHeteroNi{i} = optimizeCbModel(model,'max');
    
    if SolHeteroNi{i}.f >= 0.1 * SolHeteroWild.f
        SolHeteroNiYN{i} = SolHeteroNi{i}.f/SolHeteroWild.f ;
    else
        SolHeteroNiYN{i} = 0;
    end
    
     model = changeRxnBounds (model,['sink_',met],0,'l');
end



%% mixotrophic condition Nitrogen

model = changeRxnBounds (model, 'EX_photon(e)', -100,'l');
model = changeRxnBounds (model, 'EX_glc(e)', -1,'l');
model = changeRxnBounds (model, 'EX_hco3(e)',-3.7,'l');

for i=1:length(Carbon)
    met = Carbon{i};
    model = changeRxnBounds (model,['sink_',met],-1,'l');
    
    SolMixNi{i} = optimizeCbModel(model,'max');
    
    if SolMixNi{i}.f >= 0.1 * SolHeteroWild.f
        SolMixNiYN{i} = SolMixNi{i}.f/SolMixWild.f ;
    else
        SolMixNiYN{i} = 0;
    end
    
     model = changeRxnBounds (model,['sink_',met],0,'l');
end


