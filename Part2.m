%% Exercise 3
% Part 2

clear all

%% b) 

% Defining Parameters 

Ka = 0;
Ki = 0;
% Rates:
V1 = 0.2 ; V2 = 0.25 ; V5 = 0.75 ; V6 = 0.75 ; V9 = 0.5 ; V10 = 0.5 ;

% Michaelis-Mentem Coefficients
K_M1 = 10 ; K_M2 = 8 ; K_M3 = 15 ; K_M4 = 15 ; K_M5 = 15 ; K_M6 = 15 ; K_M7 = 15 ;
K_M8 = 15 ; K_M9 = 15 ; K_M10 = 15 ;

% Catalyze Coefficients:
k_cat3 = 0.025; k_cat4 = 0.025; k_cat7 = 0.025; k_cat8 = 0.025;

% Initial Concentrations:
MAPKKK0 = 100 ;
MAPKK0 = 300 ;
MAPK0 = 300 ;

% Initial Conditions:
X0 = zeros(8,1) ; % X = [MAPKKK ; MAPKKK_P ; MAPKK ; MAPKK_P ; MAPKK_PP ; MAPK ; MAPK_P ; MAPK_PP]
X0(1)=100;
X0(3)=300;
X0(6)=300;

% Time Span:
T_Span = [0 5000]; 


% Writing Differential Equations:

dMAPKKK = @(t, X) V2*X(2)/(K_M2 + X(2)) - V1*X(1)/(K_M1 + X(1)) ;
dMAPKKK_P = @(t, X) -V2*X(2)/(K_M2 + X(2)) + V1*X(1)/(K_M1 + X(1)) ;
dMAPKK = @(t, X) V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ;
dMAPKK_P = @(t, X) -V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ... 
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4));
dMAPKK_PP = @(t, X) -V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4));
dMAPK = @(t, X) V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6));
dMAPK_P = @(t, X) -V10*X(7)/(K_M10 + X(7)) + k_cat7*X(5)*X(6)/(K_M7 + X(6))...
    - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8));
dMAPK_PP = @(t, X) k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8));


% Defining the Function for Numerically Solve the Differential Equation
f = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)/(K_M1 + X(1)) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)/(K_M1 + X(1)); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ; 

% Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
[t1,y1] = ode45 (f, T_Span, X0) ; 


% Plotting the Results
figure 
plot(t1, y1) ;
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('Time Evolution of the Species');
hold on



%% c)

% Defining Parameters
V1max_vec = 0:0.01:0.4 ; 
Ka = 0;
Ki = 0;
MAPK_PP_concentration= [] ; 

% Solving for Different Values of V1
for i = 1:length(V1max_vec) 
    V1_prime = V1max_vec(i) ;
    
    % Defining the Function for Numerically Solve the Differential Equation
    f2 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1_prime*X(1)/(K_M1 + X(1)) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1_prime*X(1)/(K_M1 + X(1)); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ; 

    % Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
    [t2,y2] = ode45 (f2, T_Span, X0) ;
    
    MAPK_PP_concentration(i) = y2(length(t2),8) ; 
end 



% Plotting the Results
figure
plot (V1max_vec, MAPK_PP_concentration)
xlabel('V1 Values [nM s^-1]')
ylabel('Concentration of MAPK-PP (nM)')
title('Dose-Response Curve');
hold on


%% d)
%% 1: V1 = 2.5

% Defining Parameters
Ka=0;
Ki=0.1;
V1=2.5; 

T_Span = [0 5000];

% Defining the Function for Numerically Solve the Differential Equation
f3 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ;

% Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
[t3,y3] = ode45 (f3, T_Span, X0) ; 

% Plotting the Results
figure 
plot(t3,  y3) ;
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('V1 = 2.5');
hold on

%% 2: V1=1

% Defining Parameters
Ka=0;
Ki=0.1;
V1=1; 

T_Span = [0 5000];

% Defining the Function for Numerically Solve the Differential Equation
f4 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ;

% Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
[t4,y4] = ode45 (f4, T_Span, X0) ; 

% Plotting the Results
figure 
plot(t4,  y4) ;
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('V1 = 1');
hold on


%% 3: V1=0.5

% Defining Parameters
Ka=0;
Ki=0.1;
V1=0.5; 

T_Span = [0 5000];

% Defining the Function for Numerically Solve the Differential Equation
f5 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ;

% Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
[t5,y5] = ode45 (f5, T_Span, X0) ; 

% Plotting the Results
figure 
plot(t5,  y5) ;
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('V1 = 0.5');
hold on

%% 4: V1 =0.2


% Defining Parameters
Ka=0;
Ki=0.1;
V1=0.2; 

T_Span = [0 5000];

% Defining the Function for Numerically Solve the Differential Equation
f6 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ;

% Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
[t6,y6] = ode45 (f6, T_Span, X0) ; 

% Plotting the Results
figure 
plot(t6,  y6) ;
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('V1 = 0.2');
hold on

%% 5: V1 =0.1

% Defining Parameters
Ka=0;
Ki=0.1;
V1=0.1; 

T_Span = [0 5000];

% Defining the Function for Numerically Solve the Differential Equation
f7 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ;

% Using ode45 Solver in The Specified Time Span and with the Defined Initial Conditions
[t7,y7] = ode45 (f7, T_Span, X0); 

% Plotting the Results
figure 
plot(t7,  y7);
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('V1 = 0.1');
hold on


%% 6: V1 =0.25

% The oscillation stops below this value of V1.
Ka=0;
Ki=0.1;
V1=0.25; 

T_Span = [0 5000];

f8 = @(t, X) [V2*X(2)/(K_M2 + X(2)) - V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))) ; ...
    - V2*X(2)/(K_M2 + X(2)) + V1*X(1)*(1+Ka*X(8))/((K_M1 + X(1))*(1+Ki*X(8))); ...
    V6*X(4)/(K_M6 + X(4)) - k_cat3*X(2)*X(3)/(K_M3 + X(3)) ; ...
    - V6*X(4)/(K_M6 + X(4)) + k_cat3*X(2)*X(3)/(K_M3 + X(3)) + ...
    V5*X(5)/(K_M5 + X(5)) - k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    - V5*X(5)/(K_M5 + X(5)) + k_cat4*X(2)*X(4)/(K_M4 + X(4)) ; ...
    V10*X(7)/(K_M10 + X(7)) - k_cat7*X(5)*X(6)/(K_M7 + X(6)) ; ...
    - V10*X(7)/(K_M10 + X(7)) + ...
    k_cat7*X(5)*X(6)/(K_M7 + X(6)) - k_cat8*X(5)*X(7)/(K_M8 + X(7)) + V9*X(8)/(K_M9 + X(8)) ; ...
    k_cat8*X(5)*X(7)/(K_M8 + X(7)) - V9*X(8)/(K_M9 + X(8))] ;



[t8,y8] = ode45 (f8, T_Span, X0) ; 

figure 
plot(t8,  y8) ;
legend ('MAPKKK','MAPKKK-P','MAPKK','MAPKK-P','MAPKK-PP','MAPK','MAPK-P','MAPK-PP')
xlabel('Time (s)')
ylabel('Concentration (nM)')
title('V1 = 0.25');
hold on


