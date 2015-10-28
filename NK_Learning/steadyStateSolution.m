% this is steady state solution

% create cell of parameters

paramDefinition;

clear

load paramValues

eta     = paramValues.eta;
beta    = paramValues.beta;
epsilon = paramValues.epsilon;
alpha   = paramValues.alpha;
delta   = paramValues.delta;
Xst     = epsilon/(epsilon-1);
rho     = paramValues.rho;
theta   = paramValues.theta;
alphaPi = paramValues.alphaPi;
alphaY  = paramValues.alphaY;
    
% evaluate function
    
% L_init    = 1;
% Lst       = fzero(@solveLabour,L_init);
% kst       = (( (Xst/alpha)*(1/beta-1+delta) )^(1/(alpha-1)))*Lst;
% Yst       = (kst^alpha)*(Lst^(1-alpha));
% Invst     = kst*delta;
% cst       = Yst - Invst;
% rKst      = (1/beta) - 1 + delta;
% wst       = (1/Xst)*(1-alpha)*(kst^alpha)*(Lst^(1-alpha));
% Rst       = 1/beta;
% profit_st = (1-(1/Xst))*Yst;
% Tst       = cst+kst-wst*Lst-Rst*kst-profit_st; 


Lst       =  ((1-alpha)*(Xst/alpha*(1/beta-1+delta))^(alpha/(alpha-1))/(Xst*(Xst/alpha*(1/beta-1+delta)-delta)*(Xst/alpha*(1/beta-1+delta))^(1/(alpha-1))))^(1/eta);

kst       =  (Xst/alpha*(1/beta-1+delta))^(1/(alpha-1))*Lst;

Yst       =   kst^alpha*Lst^(1-alpha);

Rst       =   1/beta;

rKst      =  (1/beta) - 1 + delta;

wst       =  (1-alpha)*kst^alpha*Lst^(-alpha)/Xst;

Invst     =   kst*delta;

cst       =   Yst - Invst;

profit_st =   (1-(1/Xst))*Yst;

Tst       =   0; 



