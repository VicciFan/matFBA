function [x_new,t] = a_to_b(x,p)

% Inputs:
% x: current state ([A, B])
% p: struct of parameters

%Outputs:
% t: new time step
% x_new: new [A,B]


Na = x(1);
Nb = x(2);

% Extract Parameters
cf = p.kf;
cb = p.kb;
a1 = cf*Na;
a2 = cb*Nb;
a0 = a1 + a2;

% Define Uniform Distribution
rand1 = rand;
rand2 = rand;

% New time step
t_updated = 1/a0*log(1/rand1);

if (rand2*a0)<=a1
   Nb = Nb+1;
   Na = Na-1;
else
   Na = Na +1;
   Nb = Nb-1;
end

% Output Values 
x_new = [Na Nb];
t = t_updated;
end

