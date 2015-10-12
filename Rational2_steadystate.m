
function [ys,check]=Rational2_steadystate(ys,exe)

global M_ lgy_

global h_prim_st h_st L_st Y_st R_st m_star_st c_prim_st c_st w_st b_st b_prim_st q_st lambda_st X_st T_st F_st ratio_hh alpha_pi alpha_Y

global epsilon v gamma beta paramj eta theta m r_R r_pi r_y rho

if isfield(M_,'param_nbr') == 1
NumberOfParameters = M_.param_nbr;
for i = 1:NumberOfParameters
  paramname = deblank(M_.param_names(i,:));
  eval([ paramname ' = M_.params(' int2str(i) ');']);
end
check = 0;
end

% parameters 

epsilon=11; % price mark up parameter, retailers

v=0.03; % share of house in production function

gamma=0.98; % time preference discount parameter of impatient entrepreneur

beta=0.99; % time preference discount parameter of patient household

paramj=0.1; % parameter front house in household utility function

eta=1.01; % labour elasticity

theta=0.75; % retailer sticky price paramter, probability that he does not change his price for one period

m=0.89;  % loan to value ratio parameter

rho=0.8; % persistency parameter of TFP shock

%policy parameters, something we can try different values

r_pi=0.27; % Taylor rule parameter front inflation

r_R=0.73; % Taylor rule parameter front lag of real interest rate

r_y=0; %check it is 0 in Iacoviello's paper, Taylor rule parameter front output

alpha_pi=0.5; % macroprudential parameter,

alpha_Y=-0.5; % macroprudential parameter

% end of parameters part


% steady state solution

m_star_st=m; % steady state value of loan to value ratio, we make it to be m so that central bank do nothing with firm borrowing at steady state

ratio_hh=fzero(@find_h,0); % ratio of household housing and entrepreneur housing

h_prim_st=1/(1+ratio_hh); % household housing

h_st=1-h_prim_st; % entrepreneur housing, used for production 

L_st=((1+(1-v)*(epsilon-1))/((1-v)*(epsilon-1)-m_star_st*paramj*ratio_hh*(1-v)*(epsilon-1)))^(-1/eta); % labour

Y_st=h_st^v*L_st^(1-v); % aggregete output

R_st=1/beta; %  real interest rate

c_prim_st=(1-v)*(epsilon-1)*Y_st/epsilon/L_st^(eta); % households/lender consumption

c_st=Y_st-c_prim_st; % entrepreneur/borrower consumption

w_st=(1-v)*(epsilon-1)/epsilon*Y_st/L_st; % real wage

b_st=m_star_st*beta*paramj*c_prim_st*h_st/(1-beta)/h_prim_st; % borrowing

b_prim_st=-b_st; % lending

q_st=paramj/(1-beta)*c_prim_st/h_prim_st; % house price

lambda_st=1/c_st-gamma*R_st/c_st; % lagrangian multipiler of borrowing constraint

X_st=epsilon/(epsilon-1); % price mark up

F_st=(1-1/X_st)*Y_st;

T_st=c_prim_st+R_st*b_prim_st-b_prim_st-w_st*L_st-F_st;

%end of solutions part




%set log-linearized deviation initially equal to zero

 Y=0;
 c=0;
 c_prim=0;
 R=0;
 pi=0;
 h_prim=0;
 q=0;
 L=0;
 w=0;
 h=0;
 lambda=0;
 X=0;
 m_star=0;
 b=0;
 b_prim=0;
 T=0;
 A=0; % productivity shock

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






