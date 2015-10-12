var Y c c_prim R pi h_prim q L w h lambda X m_star b b_prim T A;

varexo e_A;

parameters h_prim_st h_st L_st Y_st R_st m_star_st c_prim_st c_st w_st b_st b_prim_st q_st lambda_st X_st T_st F_st ratio_hh epsilon v gamma beta paramj eta theta m r_pi r_y r_R chi rho alpha_pi alpha_Y;

model(linear);

//equation 1

Y_st*Y=c_st*c+c_prim_st*c_prim;

//equation 2

R-c_prim(+1)-pi(+1)=-c_prim;

//equation 3

-paramj/h_prim_st*h_prim=q_st/c_prim_st*(q-c_prim)-beta*q_st/c_prim_st*(q(+1)-c_prim(+1));

//equation 4

c_prim+(eta-1)*L=w;

//equation 6

Y=A+v*h(-1)+(1-v)*L;

A=rho*A(-1)+e_A;

//equation 7

-1/c_st*c=gamma/beta/c_st*R-gamma/beta/c_st*pi(+1)-gamma/beta/c_st*c(+1)+lambda_st*lambda;

//equation 8

q_st/c_st*(q-c)=gamma*q_st/c_st*(q(+1)-c(+1))+gamma*v*Y_st/h_st/X_st/c_st*(Y(+1)-h-X(+1)-c(+1))+lambda_st*m_star_st*q_st*beta*(lambda+m_star+q(+1)+pi(+1)-R);

//equation 9

Y_st/X_st*(Y-X)+b_st*b=c_st*c+q_st*h_st*(h-h(-1))+b_st/beta*(R(-1)+b(-1)-pi)+w_st*L_st*(w+L);

//equation 10

Y-X-L=w;

//equation 11

beta*pi(+1)=pi+(1-theta*beta)*(1-theta)/theta*X;

//equation 12

b=m_star+q(+1)+h+pi(+1)-R;

b_st*b=-b_prim_st*b_prim+q_st*h_st/R_st*m_star_st*m_star-m*q_st*h_st/R_st*(q(+1)+h+pi(+1)-R);

//equation 13

h_st*h+h_prim_st*h_prim=0;

//equation 14

T=c_prim_st*c_prim+q_st*h_prim_st*(h_prim-h_prim(-1))+R_st*b_prim_st*(R(-1)+b_prim(-1)-pi)-b_prim_st*b_prim-w_st*L_st*(w+L)-(1-1/X_st)*Y_st*Y-Y_st/X_st*X;

//equation 15

m_star=alpha_pi*pi(+1)+alpha_Y*Y;

//equation 16

R=(1-r_R)*((1+r_pi)*pi(-1)+r_y*Y(-1))+r_R*R(-1);


end;

shocks;

var e_A; // productivity shock

stderr 1;

end;


resid(1);

steady;

check;

stoch_simul(irf=0,order = 1);

conditional_variance_decomposition=1;
