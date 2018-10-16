%% Ex4

clear all;
close all;
clc;
%% Load Dataset
load('ecoli_core_model.mat');

% Choose the solver
changeCobraSolver('mosek');

% Define the biomass reaction as the objective function.
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

%% 1
%%Robustness Analysis
% Set the max import rate of O2
model = changeRxnBounds(model,'EX_o2(e)',-20,'b'); 

GrowthRate_Glucose = zeros(41,1);

% Increase Glucose uptake from 1 to 40
for i = 0:40 %number of iterations
    % Set the import boundries
    model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
    
    % Get max biomass production by optimizing the model
    OptSolution1 = optimizeCbModel(model,'max');
    
    % Find max growth rate
    GrowthRate_Glucose(i+1) = OptSolution1.f;
end

%%If Oxygen rate is constant: 
% Plot the growth yield as a function of carbon uptake
figure('Name','Glucose')
set(figure(1),'Position', [0, 900, 1300, 400]);
subplot(1,3,1)
plot(0:40, GrowthRate_Glucose)
title('Robustness Analysis Glucose')
xlabel('Glucose uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Growth rate (hr\^-1)')

%% 2 
%%Phenotypic Phase Plane
PhoneticPhases= zeros(31);
ShadowPrices_Glucose = zeros(31);

% Change Oxygen and Glucose uptake and calculate the growth rate of the objective function

for i=0:30
    for j=0:30
        % Set the import boundries
        model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
        model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
        
        % Get max biomass production by optimizing the model
        OptSolution2 = optimizeCbModel(model,'max');
        
        % Find max growth rate
        PhoneticPhases(i+1,j+1) =  OptSolution2.f;
        
        
        if OptSolution2.f > 0 % Continue updating the shadow prices
            ShadowPrices_Glucose(i+1,j+1) = OptSolution2.y(57);
        else % Set shadow price to 0 if we have no flux
            ShadowPrices_Glucose(i+1,j+1) = 0;
        end
        
   end
end

% Plot phenotypic phase plane graph
subplot(1,3,2)
surfl(PhoneticPhases.')
title('Phenotypic Phase Plane Glucose')
xlabel('Glucose uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Oxygen uptake rate (mmol gDW\^-1 hr\^-1)')
zlabel('Growth rate (hr\^-1)')

% Plot shadow price graph
subplot(1,3,3)
pcolor(ShadowPrices_Glucose.')
title('Shadow Prices for Glucose')
xlabel('Glucose uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Oxygen uptake rate (mmol gDW\^-1 hr\^-1)')

%% 3
%%Shadow Prices

% Find most repeated values among shadow prices
ShadowPrices_Glucose= round(ShadowPrices_Glucose,4);
ShadowPrices_Glucose_M = unique(ShadowPrices_Glucose);
L=length(ShadowPrices_Glucose_M);
ind= zeros(L);
ShadowPrices_Final = zeros(10,1);
for i=1:L
    % Count the number of values in ShadoPrices within each bin range
    ind(i)=histc(ShadowPrices_Glucose(:),ShadowPrices_Glucose_M(i));
    
    % If more than 30, set to max value
    if ind(i) > 30
        ShadowPrices_Final(i)= ShadowPrices_Glucose_M(i);
    end
end
ShadowPrices_Final=unique(ShadowPrices_Final);
DistincPhases_Num=length(ShadowPrices_Final);


%% 4 
%%Succinate

%%%1
model = changeRxnBounds(model,'EX_o2(e)',-20,'b');

%Reset succinate uptake rate
model = changeRxnBounds(model,'EX_pyr(e)',0,'b');
SuccinateRate = zeros(41,1);

for i = 0:40
    model = changeRxnBounds(model,'EX_succ(e)',-i,'b');
    OptSolution1 = optimizeCbModel(model,'max');
    SuccinateRate(i+1) = OptSolution1.f;
end

% Plot 
figure('Name','Succinate')
set(figure(3),'Position', [0, 1000, 1450, 400]);
subplot(1,3,1)
plot(0:40, SuccinateRate)
title('Robustness Analysis Succinate')
xlabel('Succinate uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Growth rate (hr\^-1)')

%%%2
SuccinateRate2= zeros(31);
ShadowPrices_Succinate=zeros(31);

for i=0:30
    for j=0:30
        model = changeRxnBounds(model,'EX_succ(e)',-i,'b');
        model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
        OptSolution2 = optimizeCbModel(model,'max');
        SuccinateRate2(i+1,j+1) =  OptSolution2.f;
        if OptSolution2.f > 0
            ShadowPrices_Succinate(i+1,j+1) = OptSolution2.y(57);
        else
            ShadowPrices_Succinate(i+1,j+1) = 0;
        end
    end
end

% Plot
subplot(1,3,2) 
surfl(SuccinateRate2.')
title('Phenotypic Phase Plane Succinate')
xlabel('Succinate uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Oxygen uptake rate (mmol gDW\^-1 hr\^-1)')
zlabel('Growth rate (hr\^-1)')
subplot(1,3,3)
pcolor(ShadowPrices_Succinate.')
title('Shadow Prices for Succinate')
xlabel('Succinate uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Oxygen uptake rate (mmol gDW\^-1 hr\^-1)')

%%%3
ShadowPrices_Succinate= round(ShadowPrices_Succinate,4);
ShadowPrices_Succinate_M = unique(ShadowPrices_Succinate);
L=length(ShadowPrices_Succinate_M);
ind= zeros(L);
ShadowPrices_SuccinateFinal = zeros(10,1);
for i=1:L
    ind(i)=histc(ShadowPrices_Succinate(:),ShadowPrices_Succinate_M(i));
    if ind(i) > 80 %Max boundry
        ShadowPrices_SuccinateFinal(i)= ShadowPrices_Succinate_M(i);
    end
end
ShadowPrices_SuccinateFinal=unique(ShadowPrices_SuccinateFinal);
DistincPhases_Num_Succinate=length(ShadowPrices_SuccinateFinal);


%%Pyruvate

%%%1
model = changeRxnBounds(model,'EX_o2(e)',-20,'b');

%Reset Pyruvate uptake rate
model = changeRxnBounds(model,'EX_pyr(e)',0,'b');
PyruvateRate = zeros(41,1);

for i = 0:40
    model = changeRxnBounds(model,'EX_succ(e)',-i,'b');
    OptSolution1 = optimizeCbModel(model,'max');
    PyruvateRate(i+1) = OptSolution1.f;
end

% Plot 
figure('Name','Pyruvate')
set(figure(3),'Position', [0, 1000, 1450, 400]);
subplot(1,3,1)
plot(0:40, PyruvateRate)
title('Robustness Analysis Pyruvate')
xlabel('Pyruvate uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Growth rate (hr\^-1)')

%%%2
PyruvateRate2= zeros(31);
ShadowPrices_Pyruvate=zeros(31);

for i=0:30
    for j=0:30
        model = changeRxnBounds(model,'EX_succ(e)',-i,'b');
        model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
        OptSolution2 = optimizeCbModel(model,'max');
        PyruvateRate2(i+1,j+1) =  OptSolution2.f;
        if OptSolution2.f > 0
            ShadowPrices_Pyruvate(i+1,j+1) = OptSolution2.y(57);
        else
            ShadowPrices_Pyruvate(i+1,j+1) = 0;
        end
    end
end

% Plot
subplot(1,3,2) 
surfl(PyruvateRate2.')
title('Phenotypic Phase Plane Pyruvate')
xlabel('Pyruvate uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Oxygen uptake rate (mmol gDW\^-1 hr\^-1)')
zlabel('Growth rate (hr\^-1)')
subplot(1,3,3)
pcolor(ShadowPrices_Pyruvate.')
title('Shadow Prices for Pyruvate')
xlabel('Pyruvate uptake rate (mmol gDW\^-1 hr\^-1)')
ylabel('Oxygen uptake rate (mmol gDW\^-1 hr\^-1)')

%%%3
ShadowPrices_Pyruvate= round(ShadowPrices_Pyruvate,4);
ShadowPrices_Pyruvate_M = unique(ShadowPrices_Pyruvate);
L=length(ShadowPrices_Pyruvate_M);
ind= zeros(L);
ShadowPrices_PyruvateFinal = zeros(10,1);
for i=1:L
    ind(i)=histc(ShadowPrices_Pyruvate(:),ShadowPrices_Pyruvate_M(i));
    if ind(i) > 80 %Max boundry
        ShadowPrices_PyruvateFinal(i)= ShadowPrices_Pyruvate_M(i);
    end
end
ShadowPrices_PyruvateFinal=unique(ShadowPrices_PyruvateFinal);
DistincPhases_Num_Pyruvate=length(ShadowPrices_PyruvateFinal);


