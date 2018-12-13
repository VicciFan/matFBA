%% Exercise 2
% Part 1

% Determining Constants
K_D = 0.25;
K_I = 1;
V_max = 50;

S = linspace(0,10,100);
v = zeros(100,1);

% Different [I] Values
I = [0 1 2 3];

figure
for k=I(1): I(4)

        f = @(s) V_max/(1+(k*K_D/s*K_I)+(k/K_I)+(K_D/s));
        syms s;
        Km = solve(V_max/(1+(k*K_D/s*K_I)+(k/K_I)+(K_D/s)) == (f(100)/2), s);
        
        %v(i) = V_max/(K_D*(S(i)+1/K_D+k/(S(i)*K_I)+k/(K_D*K_I)));
        
    fplot(f,[0,10]);
    hold on
end

xlabel('Substrate concentration [mM]')
ylabel('Reaction rate [mmol/min]')
title('Sunctrate Concentraions');
legend('[I] = 0 mM','[I] = 1 mM','[I] = 2 mM','[I] = 3 mM')
hold on

