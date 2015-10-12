% hello

% this is the main file for households learning in Iacoviello (2015), with dynare and all the other stuff

% first thing. Run dynare file and obtain transition matrix

dynare Rational2.mod % This is dynare file with RE 

[transitionVar_A, transitionShock_B] = transition_matrix(oo_); 
% this function takes the transition matrix A and B from the dynare result,
% A is for variables and B is for shocks
transitionVar_A(7,:)=transitionVar_A(8,:);

transitionVar_A=transitionVar_A(1:7,:); 
% We have this two line syntax because this function does not work perfectly
% one row are all zero, so we need to remove it from matrix A
clearvars -except transitionVar_A transitionShock_B 
% clear all the variables except the transition matrixes 

time_period=40; % time period of learning

forecast_period=100; 
% time period that households make forecasts. The larger it is, the more precise the result of consumption under PLM will be

% constant gain learning section

gainParam = 0.04; % gain parameter, gamma in tex file.

% initialize the predetermined variables at the steady states to be equal to zero at time 0 (first column). 
%(log-deviation, except the Central bank money transfer T is the absolute deviation)

All_variables=zeros(length(transitionShock_B),time_period); 
% space for all variables, first column is for time 0, second for time 1, so on so forth.

Pred_variables=zeros(size(transitionVar_A,1), time_period);
% space for predetermined variables, first column is for time 0, second for time 1, so on so forth.
% The predeterminded variables are interest rate R (4th), households's housing h_prim (6th), entrepreneur's housing h (10th),  borrowing b (14th), lending b_prim (15th), Technological shock A (17th), inflation pi (5th).

% initialize the shock

e_A=0.01*randn; % shock happens from time 0 to 1.

rho=0.8; % decaying rate of shock

Steadystate_file;
% get the steady state solutions, call 'Steadystate_file.m'

% Below shows consumption equation parameters, details of parameters are in tex file in Appendix (7.4 Households learning structure).
% The notation of the parameters are the same as in tex file.

C_1=c_prim_st*(1+paramj*(1+beta)/(1-beta)^2+1/(eta-1)*w_st^(eta/(eta-1))*c_prim_st^(-eta/(eta-1))+(1-v)/(eta-1)*(1-1/X_st)*(1-paramj*c_prim_st/(1-beta)/q_st)^v*(w_st/c_prim_st)^((1-v)/(eta-1)-1)*w_st/c_prim_st^2-v*paramj*beta*(1-1/X_st)*(w_st/c_prim_st)^((1-v)/(eta-1))*(1-paramj*c_prim_st/(1-beta)/q_st)^(v-1)/q_st/(1-beta)^2);

C_2=-paramj*beta/(1-beta)^2*c_prim_st;

C_0=c_prim_st*(-paramj/(1-beta)^2+v*paramj*(1-1/X_st)*(w_st/c_prim_st)^((1-v)/(eta-1))*(1-paramj*c_prim_st/q_st/(1-beta))^(v-1)/q_st/(1-beta)^2);

C_q1=q_st*(-paramj*c_prim_st*(1+beta)/(1-beta)^2/q_st+v*paramj*beta*(1-1/X_st)*(w_st/c_prim_st)^((1-v)/(eta-1))*(1-paramj*c_prim_st/q_st/(1-beta))^(v-1)/q_st^2/(1-beta)^2*c_prim_st);

C_q2=paramj*beta*c_prim_st/(1-beta)^2;

C_q0=q_st*(paramj*c_prim_st/(1-beta)^2/q_st-v*paramj*(1-1/X_st)*(w_st/c_prim_st)^((1-v)/(eta-1))*(1-paramj*c_prim_st/(1-beta)/q_st)^(v-1)*c_prim_st/(1-beta)^2/q_st^2);

C_X0=-(1-paramj*c_prim_st/(1-beta)/q_st)^v*(w_st/c_prim_st)^((1-v)/(eta-1))/X_st;

C_w=-w_st*(eta/(eta-1)*c_prim_st^(-1/(eta-1))*w_st^(1/(eta-1))+(1-v)/(eta-1)*(1-1/X_st)*(1-paramj*c_prim_st/(1-beta)/q_st)^v*(w_st/c_prim_st)^((1-v)/(eta-1)-1)/c_prim_st);

H_pi(1)=chi_st*(1/R_st)/(1-1/R_st)-theta/(1-theta*beta)/(1-theta)*C_X0/R_st;

for i=2:forecast_period

H_pi(i)=chi_st*(1/R_st)^i/(1-1/R_st);

end 

H_R0=chi_st*(1/R_st)/(1-1/R_st);

for i=1:forecast_period

H_R(i)=chi_st*(1/R_st)^(i+1)/(1-1/R_st);

end

for i=1:forecast_period
    
   H_w(i)=C_w*(1/R_st)^i;
    
end

for i=1:forecast_period
    
   H_T(i)=(1/R_st)^i;
    
end

for i=1:forecast_period
    
   H_X(i)=C_X0*(1/R_st)^i;
    
end

H_c0=1/R_st*C_0+C_1;

H_c(1)=(1/R_st)^2*C_0+1/R_st*C_1+C_2;

for i=2:forecast_period
    
   H_c(i)=(1/R_st)^(i+1)*C_0+(1/R_st)^i*C_1+(1/R_st)^(i-1)*C_2;
    
end

H_q0=1/R_st*C_q0+C_q1;

H_q(1)=(1/R_st)^2*C_q0+(1/R_st)*C_q1+C_q2;

for i=2:forecast_period

H_q(i)=(1/R_st)^(i+1)*C_q0+(1/R_st)^i*C_q1+(1/R_st)^(i-1)*C_q2;

end

% The end of learning parameters.

 
phi_b_prim(:,1)  = 0.001*randn(3,1); %initialize initial parameters of PLM in constant gain learning, phi_(variable names) are parameters estimated in PLM

phi_Infl(:,1) = 0.001*randn(3,1); 
               
phi_Intr(:,1) =  0.001*randn(3,1);
               
phi_q(:,1) =0.001*randn(3,1); 
               
phi_w(:,1) =0.001*randn(3,1);

phi_T(:,1) =0.001*randn(3,1);

D_Old =  eye(3);  % initialize D matrix in constant gain learning (5.1 PLM of households)

% where loop starts

for t=2:time_period

if t==2
    
    All_variables(:,t)=transitionVar_A'*Pred_variables(:,t-1)+transitionShock_B'*e_A;
    % from t-1 to t, the shock comes
else 
    
    All_variables(:,t)=transitionVar_A'*Pred_variables(:,t-1); 
    % Technological parameter is included in All_variables, so we just shock it one time and A has decaying rate rho
    
end
% new information is needed to update parameter and D matrix in PLM,
% we get them from RE matrix at time 1.

w_new=(1+All_variables(9,t))*w_st; % wage is 9th variable

pi_new=1+All_variables(5,t); % inflation is 5th variable

q_new=q_st*(1+All_variables(7,t));% house price is 7th variable

R_new=(All_variables(4,t)+1)*R_st; % interest rate is 4th variable

T_new=All_variables(16,t);

% we take the information needed for learning: interest rate, inflation,
% wage rate and house prices, into the learning prodecure

% people choose parameter values arbitrately at the beginning

zMat = [1 All_variables(15,t-1) All_variables(17,t-1) ]'; % initialize state variables z=(1 b_prim_hat A_hat) at time t-1 with intercept


% CG Learning Algorithm, phi_(variable name) are parameters estimated. Equations are at section 'PLM of households,prediction' 

    D_New = D_Old + gainParam * (zMat * zMat' - D_Old);    % update D matrix
    
    % update parameters, Equation in Section 'PLM of households'
    
    phi_b_prim(:,t)    = phi_b_prim(:,t-1) + gainParam * D_New \ zMat * (All_variables(15,t) - phi_b_prim(:,t-1)' * zMat);
    
    phi_Infl(:,t) = phi_Infl(:,t-1) + gainParam * D_New \ zMat * (All_variables(5,t) - phi_Infl(:,t-1)' * zMat);
    
    phi_Intr(:,t) = phi_Intr(:,t-1) + gainParam * D_New \ zMat * (All_variables(4,t) - phi_Intr(:,t-1)' * zMat);
    
    phi_q(:,t)    = phi_q(:,t-1) + gainParam * D_New\ zMat * ( All_variables(7,t)- phi_q(:,t-1)' * zMat);
    
    phi_w(:,t)    = phi_w(:,t-1) + gainParam * D_New\ zMat * (All_variables(9,t) - phi_w(:,t-1)' * zMat);
    
    phi_T(:,t)    = phi_w(:,t-1) + gainParam * D_New\ zMat * (All_variables(16,t) - phi_T(:,t-1)' * zMat);
    
    D_Old=D_New; 
    
    zMat=[1 All_variables(15,t) All_variables(17,t)]'; % update regressors information
   
    
    % Below is prediction of variables after households update their
    % parameters. Equation 'Prediction' in tex file.
    
    Prediction_b_prim_hat=zeros(3,forecast_period); % we have [1 b_prim, A_hat] each column, first column is 1 step ahead forcast, so on so forth
    
    Prediction_PLM_hat=zeros(5,forecast_period);  % prediction of inflation, interest rate, house prices, wage and public transfers
    
    Prediction_b_prim_hat(1,1)=1; % intercept
    
    Prediction_b_prim_hat(2,1)=phi_b_prim(:,t)'*zMat; % 1 step ahead forecast of b_prim, we put it in first column
    
    Prediction_b_prim_hat(3,1)=rho*All_variables(17,t);
    
    Prediction_PLM_hat(1,1)= phi_T(:,t)'*zMat;
    
    Prediction_PLM_hat(2,1)=phi_Infl(:,t)'*zMat; % second row of prediction matrix is inflation 
   
    Prediction_PLM_hat(3,1)=phi_Intr(:,t)'*zMat; % third row of prediction matrix is interest rate 
   
    Prediction_PLM_hat(4,1)=phi_q(:,t)'*zMat; % fourth row of prediction matrix is house price
   
    Prediction_PLM_hat(5,1)=phi_w(:,t)'*zMat; % fifth row of prediction matrix is wage

   
   % make n-step ahead forecast, used in households choosing consumption.
    
for i=2:forecast_period;
    
    Prediction_b_prim_hat(1,i)=1;
    
    Prediction_b_prim_hat(2,i)=phi_b_prim(:,t)'*Prediction_b_prim_hat(:,i-1);
    
    Prediction_b_prim_hat(3,i)=rho*Prediction_b_prim_hat(2,i-1);
    
   Prediction_PLM_hat(1,i)=phi_T(:,t)'*Prediction_b_prim_hat(:,i-1);  
    
   Prediction_PLM_hat(2,i)=phi_Infl(:,t)'*Prediction_b_prim_hat(:,i-1);
   
   Prediction_PLM_hat(3,i)=phi_Intr(:,t)'*Prediction_b_prim_hat(:,i-1);
   
   Prediction_PLM_hat(4,i)=phi_q(:,t)'*Prediction_b_prim_hat(:,i-1);
   
   Prediction_PLM_hat(5,i)=phi_w(:,t)'*Prediction_b_prim_hat(:,i-1);
 
   
  
end 

% this calculate the double sum values in consumption function at the end
% of section 'households learning', one equation before 'Prediction'.

sumR=0; sumpi=0;

for i=2:forecast_period

sumR= sumR+H_c0*All_variables(4,t)+H_c(1:i-1)*Prediction_PLM_hat(3,1:i-1)';
    
sumpi=sumpi+H_c(1:i)*Prediction_PLM_hat(2,1:i)';

end

% this is consumption function to get current period consumption,
% log-deviation
% The equation just before section 'PLM of households'

c_prim_hat=1/(-H_c0-sum(H_c))*(R_st*b_prim_st*All_variables(5,t)+H_pi*Prediction_PLM_hat(2,:)'-H_R0*All_variables(4,t)-H_R*Prediction_PLM_hat(3,:)'+C_w*All_variables(9,t)+H_w*Prediction_PLM_hat(5,:)'-All_variables(16,t)-H_T*Prediction_PLM_hat(1,:)'+H_c(1)*(All_variables(4,t)-Prediction_PLM_hat(2,1))+sumR-sumpi+H_q0*All_variables(7,t)+H_q(1)*Prediction_PLM_hat(4,1)+H_q(2:end)*Prediction_PLM_hat(4,2:end)');

c_prim_new=c_prim_st*(1+c_prim_hat); % consumption decisions, 

L_prim_new=(w_new/c_prim_new)^(1/(eta-1)); 
% labour supply, given consumption and wage rate. 

 % we get expectation of consumption tomorrow from euler equation, so as to get the households' house holding

c_prim1_hat=c_prim_hat+All_variables(4,t)-All_variables(5,t);

h_prim_new=paramj/((All_variables(7,t-1)+1)/c_prim_new-beta*(Prediction_PLM_hat(4,1)+1)*q_st/(c_prim1_hat+1)/c_prim_st);

h_prim_new_hat=(h_prim_new-h_prim_st)/h_prim_st; 
% household's housing, log-deviation

h_new_hat=-h_prim_new_hat;
% entrepreneur's housing, log-deviation

Y_new=All_variables(17,t)*h_prim_new^v*L_prim_new^(1-v);

X_new=(1-v)*Y_new/w_new/L_prim_new;

F_new=(1-1/X_new)*Y_new;

b_prim_new=c_prim_new+q_new*(h_prim_new-(All_variables(15,t-1)+1)*h_prim_st)+(1+All_variables(4,t-1))*R_st*(1+All_variables(15,t-1))*b_prim_st/(1+All_variables(5,t))-w_new*L_prim_new-F_new-All_variables(16,t); % we derive the lending of household from intertemporal constraint, the last term is transfer from RE

% revise all the other variables in predetermined variable vector, before
% we move on to next period. 
% The predeterminded variables are interest rate R (4th), households's housing h_prim (6th), entrepreneur's housing h (10th),  borrowing b (14th), lending b_prim (15th), Technological shock A (17th), inflation pi (5th).


Pred_variables(1,t)= All_variables(4,t);% interest rate

Pred_variables(2,t) = h_prim_new_hat;

Pred_variables(3,t)=h_new_hat;

Pred_variables(5,t)=(b_prim_new-b_prim_st)/b_prim_st; % households lending

Pred_variables(4,t)=-Pred_variables(5,t)% firm's borrowing

Pred_variables(6,t)=All_variables(17,t);

Pred_variables(7,t)=All_variables(5,t);

end










