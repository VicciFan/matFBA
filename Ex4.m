% Exercise 4

close all
clear all
format long

%% a
% The function is implemented and used in this script

%% b & c
% Define final Results Matrix for All the Runs
Res = zeros(10,3,10000);

% Set Parameters
p.kf = 2;
p.kb = 1;

Volume = 0.17e-11 ;
A_Mol = 1e-9;
B_Mol = 0e-9;

% Start time
p.t = 0;    

% Initial Number of Molecules
Mol_No = 6.022140857e23;

% Simulate the Reaction for 10 Times
figure
for i=1:10
   
    A_Val = A_Mol*Mol_No*Volume;
    B_Val = B_Mol*Mol_No*Volume;
    
    % Initial Conditions
    Init_Cond = [A_Val B_Val];
    
    % Plot the Results in a Step Plot
    X = zeros(2,10000);
    T = zeros(1,10000);
    
    % Simulate for 10000 time steps
    t_new=0;
    T_tot=0;
    
    for j=1:10000       
        % Use a_to_b Function
        [Init_Cond,t]=a_to_b(Init_Cond,p);        
        % Update Values
        t_new = t_new+t;
        T_tot = T_tot+t;
        X(:,j) = Init_Cond(:);
        T(1,j) = t_new;
    end
    
    Tmp(1,:) = T(:);
    Tmp(2:3,:) = X(:,:);
    Res(i,:,:)= Tmp(:,:);
    
    % Plot the Results
    stairs(T,X(2,:));
    hold on
    stairs(T,X(1,:));
    hold on
end
xlabel('Time (s)')
ylabel('Molecules Number')
title('Simulated Reactions:Volume = 1.7pL')
legend('A', 'B') 

%% d
% Calculate Mean and Std
Tmp = Res(:,1,10000);
T_Max = min(Tmp);
T_Max = round(T_Max);

% Choose the Time Grid 
No_Cell = 100;
Step_Size = T_Max/No_Cell;
Grid = 0:Step_Size:T_Max;

% Performing the Simulations Using the Specified Grid
Tmp = zeros(3,length(Grid));
Res2 = zeros(10,3,length(Grid));
New_Data = zeros(2,length(Grid));
data = zeros(3,10000);

% Plot the Stochastic Simulations
figure
for i=1:10
    data(:,:) = Res(i,:,:);
    
    % Interpolation
    New_Data(1,:) = interp1(data(1,:), data(2,:),Grid, 'previous');
    New_Data(2,:) = interp1(data(1,:), data(3,:),Grid, 'previous');
    
    Tmp(1,:) = Grid(:);
    Tmp(2:3,:) = New_Data(1:2,:);
    
    Res2(i,:,:) = Tmp(:,:);
    
    % Plot the Results on the Step Plot
    
    stairs(Grid,New_Data(2,:));
    hold on
    stairs(Grid,New_Data(1,:));
    hold on
end

xlabel('Time (s)')
ylabel('Molecules Number')
title('Mean & Std: Volume = 1.7pL')
legend('A', 'B') 

% Calculate Mean and Std for Points on the Grid 
Data_Mean = zeros(3,length(Grid));
Data_Std = zeros(2,length(Grid));
Data_Mean(1,:) = Grid(:);
Data_Mean(2,:) = mean(Res2(:,2,:));
Data_Mean(3,:) = mean(Res2(:,3,:));
Data_Std(1,:) = std(Res2(:,2,:));
Data_Std(2,:) = std(Res2(:,3,:));

% Plot Mean and Std Curves
figure
errorbar(Grid,Data_Mean(2,:),Data_Std(1,:))
hold on
errorbar(Grid,Data_Mean(3,:),Data_Std(2,:))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Mean & Std of Grid Points: Volume = 1.7pL')
legend('A', 'B') 

%% e
% Create ODE Simulation

% Define Initial Concentrations
X_I = zeros(2,1);
X_I(1) = 1e-9; 
X_I(2) = 0e-9;

% ODE
T_Span = [0 T_tot];
f = @(time,X)[-p.kf*X(1)+p.kb*X(2);p.kf*X(1)- p.kb*X(2)] ;
opts = odeset('AbsTol',1e-12);
[time,Y] = ode45(f,T_Span, X_I, opts);
figure
plot(time,(Y*Mol_No*Volume)) 
xlabel('Time (s)')
ylabel('Molecules Number')
title('Simulated Reactions Using ODE: Volume = 1.7pL')
legend('A', 'B') 

% Plot
figure
errorbar(Grid,Data_Mean(2,:),Data_Std(1,:))
hold on
errorbar(Grid,Data_Mean(3,:),Data_Std(2,:))
hold on
plot(time,(Y*Mol_No*Volume))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Mean & Std with Results of ODE: Volume = 1.7pL')
legend('A', 'B') 

% Compare Deterministic and Stochastic
% For A
figure
errorbar(Grid,Data_Mean(2,:),Data_Std(1,:))
hold on
plot(time,(Y(:,1)*Mol_No*Volume))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Comparison for A reactant: Vol =1.7pL')



% For B
figure
errorbar(Grid,Data_Mean(3,:),Data_Std(2,:))
hold on
plot(time,(Y(:,2)*Mol_No*Volume))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Comparison for B reactant: Vol =1.7pL')


%% f
% Repeat for 0.17 pL
clear all

% Define final Results Matrix for All the Runs
Res = zeros(10,3,10000);

% Set Parameters
p.kf = 2;
p.kb = 1;

Volume = 0.17e-12 ;
A_Mol = 1e-9;
B_Mol = 0e-9;

% Start time
p.t = 0;    

% Initial Number of Molecules
Mol_No = 6.022140857e23;

% Simulate the Reaction for 10 Times
figure
for i=1:10
   
    A_Val = A_Mol*Mol_No*Volume;
    B_Val = B_Mol*Mol_No*Volume;
    
    % Initial Conditions
    Init_Cond = [A_Val B_Val];
    
    % Plot the Results in a Step Plot
    X = zeros(2,10000);
    T = zeros(1,10000);
    
    % Simulate for 10000 time steps
    t_new=0;
    T_tot=0;
    
    for j=1:10000       
        % Use a_to_b Function
        [Init_Cond,t]=a_to_b(Init_Cond,p);        
        % Update Values
        t_new = t_new+t;
        T_tot = T_tot+t;
        X(:,j) = Init_Cond(:);
        T(1,j) = t_new;
    end
    
    Tmp(1,:) = T(:);
    Tmp(2:3,:) = X(:,:);
    Res(i,:,:)= Tmp(:,:);
    
    % Plot the Results
    stairs(T,X(2,:));
    hold on
    stairs(T,X(1,:));
    hold on
end
xlabel('Time (s)')
ylabel('Molecules Number')
title('Simulated Reactions:Volume = 0.17pL')
legend('A', 'B') 

% Calculate Mean and Std
Tmp = Res(:,1,10000);
T_Max = min(Tmp);
T_Max = round(T_Max);

% Choose the Time Grid 
No_Cell = 100;
Step_Size = T_Max/No_Cell;
Grid = 0:Step_Size:T_Max;

% Performing the Simulations Using the Specified Grid
Tmp = zeros(3,length(Grid));
Res2 = zeros(10,3,length(Grid));
New_Data = zeros(2,length(Grid));
data = zeros(3,10000);

% Plot the Stochastic Simulations
figure
for i=1:10
    data(:,:) = Res(i,:,:);
    
    % Interpolation
    New_Data(1,:) = interp1(data(1,:), data(2,:),Grid, 'previous');
    New_Data(2,:) = interp1(data(1,:), data(3,:),Grid, 'previous');
    
    Tmp(1,:) = Grid(:);
    Tmp(2:3,:) = New_Data(1:2,:);
    
    Res2(i,:,:) = Tmp(:,:);
    
    % Plot the Results on the Step Plot
    
    stairs(Grid,New_Data(2,:));
    hold on
    stairs(Grid,New_Data(1,:));
    hold on
end

xlabel('Time (s)')
ylabel('Molecules Number')
title('Mean & Std: Volume = 0.17pL')
legend('A', 'B') 

% Calculate Mean and Std for Points on the Grid 
Data_Mean = zeros(3,length(Grid));
Data_Std = zeros(2,length(Grid));
Data_Mean(1,:) = Grid(:);
Data_Mean(2,:) = mean(Res2(:,2,:));
Data_Mean(3,:) = mean(Res2(:,3,:));
Data_Std(1,:) = std(Res2(:,2,:));
Data_Std(2,:) = std(Res2(:,3,:));

% Plot Mean and Std Curves
figure
errorbar(Grid,Data_Mean(2,:),Data_Std(1,:))
hold on
errorbar(Grid,Data_Mean(3,:),Data_Std(2,:))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Mean & Std of Grid Points: Volume = 0.17pL')
legend('A', 'B') 

% Create ODE Simulation

% Define Initial Concentrations
X_I = zeros(2,1);
X_I(1) = 1e-9; 
X_I(2) = 0e-9;

% ODE
T_Span = [0 T_tot];
f = @(time,X)[-p.kf*X(1)+p.kb*X(2);p.kf*X(1)- p.kb*X(2)] ;
opts = odeset('AbsTol',1e-12);
[time,Y] = ode45(f,T_Span, X_I, opts);
figure
plot(time,(Y*Mol_No*Volume)) 
xlabel('Time (s)')
ylabel('Molecules Number')
title('Simulated Reactions Using ODE: Volume = 0.17pL')
legend('A', 'B') 

% Plot
figure
errorbar(Grid,Data_Mean(2,:),Data_Std(1,:))
hold on
errorbar(Grid,Data_Mean(3,:),Data_Std(2,:))
hold on
plot(time,(Y*Mol_No*Volume))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Mean & Std with Results of ODE: Volume = 0.17pL')
legend('A', 'B') 

% Compare Deterministic and Stochastic
% For A
figure
errorbar(Grid,Data_Mean(2,:),Data_Std(1,:))
hold on
plot(time,(Y(:,1)*Mol_No*Volume))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Comparison for A reactant: Vol =0.17pL')



% For B
figure
errorbar(Grid,Data_Mean(3,:),Data_Std(2,:))
hold on
plot(time,(Y(:,2)*Mol_No*Volume))
xlabel('Time (s)')
ylabel('Molecules Number')
title('Comparison for B reactant: Vol =0.17pL')







