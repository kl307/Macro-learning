% This file solves the stteady state value of Iacoviello 2005

%Parameters

global epsilon v gamma beta paramj eta theta m r_pi m_star_st

epsilon=11;

v=0.03;

gamma=0.98;

beta=0.99;

paramj=0.1;

eta=1.01;

theta=0.75;

m=0.89;

r_pi=0.27;

m_star_st=m;

%find the steady state value of h

ratio_hh=fzero(@find_h,0);

h_prim_st=1/(1+ratio_hh);

h_st=1-h_prim_st;

L_st=((1+(1-v)*(epsilon-1))/((1-v)*(epsilon-1)-m_star_st*paramj*ratio_hh*(1-v)*(epsilon-1)))^(-1/eta);

Y_st=h_st^v*L_st^(1-v);

R_st=1/beta;

c_prim_st=(1-v)*(epsilon-1)*Y_st/epsilon/L_st^(eta);

c_st=Y_st-c_prim_st;

w_st=(1-v)*(epsilon-1)/epsilon*Y_st/L_st;

b_st=m_star_st*beta*paramj*c_prim_st*h_st/(1-beta)/h_prim_st;

b_prim_st=-b_st;

q_st=paramj/(1-beta)*c_prim_st/h_prim_st;

lambda_st=1/c_st-gamma*R_st/c_st;

X_st=epsilon/(epsilon-1);

F_st=(1-1/X_st)*Y_st;

T_st=c_prim_st+R_st*b_prim_st-b_prim_st-w_st*L_st-F_st;

chi_st=c_prim_st-w_st*L_st-F_st-T_st;

