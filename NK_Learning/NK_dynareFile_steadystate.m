
function [ys,check]=NK_dynareFile_steadystate(ys,exe)

global M_ lgy_

global Lst kst Yst Invst cst rKst wst Rst profit_st

global eta beta epsilon alpha delta Xst theta rho

if isfield(M_,'param_nbr') == 1
NumberOfParameters = M_.param_nbr;
for i = 1:NumberOfParameters
  paramname = deblank(M_.param_names(i,:));
  eval([ paramname ' = M_.params(' int2str(i) ');']);
end
check = 0;
end

% parameters 

load paramValues

eta     = paramValues.eta;
beta    = paramValues.beta;
epsilon = paramValues.epsilon;
alpha   = paramValues.alpha;
delta   = paramValues.delta;
Xst     = paramValues.Xst;
theta   = paramValues.theta;
rho     = paramValues.rho;

clear paramValues
% end of parameters part

% steady state solution
% steady state solution requires 'solveLabour.m' function

Lst       =  ((1-alpha)*(Xst/alpha*(1/beta-1+delta))^(alpha/(alpha-1))/(Xst*(Xst/alpha*(1/beta-1+delta)-delta)*(Xst/alpha*(1/beta-1+delta))^(1/(alpha-1))))^(1/eta);

kst       =  (Xst/alpha*(1/beta-1+delta))^(1/(alpha-1))*Lst;

Yst       =   kst^alpha*Lst^(1-alpha);

Rst       =   1/beta;

rKst      =  (1/beta) - 1 + delta;

wst       =  (1-alpha)*kst^alpha*Lst^(-alpha)/Xst;

Invst     =   kst*delta;

cst       =   Yst - Invst;

profit_st =   (1-(1/Xst))*Yst;



%end of solutions part

%set log-linearized deviation initially equal to zero
 
R      = 0;
w      = 0;
c      = 0;
L      = 0;
rK     = 0;
k      = 0;
Y      = 0;
Inv    = 0;
X      = 0;
infl   = 0;
A      = 0; % productivity shock
T      = 0;
profit = 0;

% This steady state variables have been checked match the dynare file.



for iter = 1:length(M_.params)
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

if isfield(M_,'param_nbr') == 1

if isfield(M_,'orig_endo_nbr') == 1
NumberOfEndogenousVariables = M_.orig_endo_nbr;
else
NumberOfEndogenousVariables = M_.endo_nbr;
end
ys = zeros(NumberOfEndogenousVariables,1);
for i = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(i,:));
  eval(['ys(' int2str(i) ') = ' varname ';']);
end
else
ys=zeros(length(lgy_),1);
for i = 1:length(lgy_)
    ys(i) = eval(lgy_(i,:));
end
check = 0;
end
end






