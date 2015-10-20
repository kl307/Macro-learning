% this is learning code for NK Model

% 1. Run this file to obtain steady state solution and parameter values

steadyStateSolution;

clear paramValues Linit

timePeriod      = 40;   % time period of model
forecastPeriod  = 50;  % forecasting period
gainParam       = 0.04; % constant gain parameter
numVars         = 12;   % number of variables in the model, without shock
allVariables    = zeros(numVars,timePeriod); % all variables at time T
CGpredictedW    = zeros(1,timePeriod); % these are predicted values from the learning algorithm
CGpredictedR    = zeros(1,timePeriod);
CGpredictedPi   = zeros(1,timePeriod);
CGpredictedT    = zeros(1,timePeriod);
CGpredictedK    = zeros(1,timePeriod);

% create vector of shocks

A_hat      = zeros(1,timePeriod+forecastPeriod);
A_hat(1,1) = 1; % initialize A_hat at 1

for j = 2:timePeriod+forecastPeriod
   
    A_hat(1,j) = A_hat(1,j-1) * 0.1*randn;
    
end

clear j

% create vector of capital

Capital      = zeros(1,timePeriod+forecastPeriod);
Capital(1,1) = kst; % initialize capital at steady state level

% how does the capital evolve then?

for j = 2:timePeriod+forecastPeriod
   
    A_hat(1,j) = A_hat(1,j-1) * 0.1*randn;
    
end

clear j


% define some parameters for consumption function

C_w       = wst*( (eta/(eta-1))*(wst^(1/(eta-1)))*(cst^(-1/(eta-1)))+((1-alpha)/(eta-1))*(1-(1/Xst))*(kst^alpha)*(cst^(-(1-alpha)/(eta-1)))*(wst^((1-alpha)/(eta-1)-1)) );
C_0       = cst*( -(eta/(eta-1))*(wst^(eta/(eta-1)))*(cst^(-eta/(eta-1)))-((1-alpha)/(eta-1))*(1-(1/Xst))*(kst^alpha)*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1)-1)) -1 );
C_X       = (1/Xst)*( (kst^alpha)*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1))) );
C_A       = (1-(1/Xst))*(kst^alpha)*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1)));
C_k       = alpha*(1-(1/Xst))*(kst^(1-alpha))*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1)));
BigLambda = wst*Lst+Tst+profit_st-cst;

% calculate sums now

G_C_temp = zeros(1,forecastPeriod);
G_w_temp = zeros(1,forecastPeriod);
G_X_temp = zeros(1,forecastPeriod);
G_A_temp = zeros(1,forecastPeriod);
G_k_temp = zeros(1,forecastPeriod);
G_T_temp = zeros(1,forecastPeriod);

% fill in the sums now

for j = 1:forecastPeriod
    
   G_C_temp(1,j) = ((1/Rst)^j)*C_0;
   G_w_temp(1,j) = ((1/Rst)^j)*C_w;
   G_X_temp(1,j) = ((1/Rst)^j)*C_X;
   G_A_temp(1,j) = ((1/Rst)^j)*C_A;
   G_k_temp(1,j) = ((1/Rst)^j)*C_k;
   G_T_temp(1,j) = ((1/Rst)^j)*Tst;
   
end

clear j

% calculate final sums

G_C = zeros(1,forecastPeriod);
G_w = zeros(1,forecastPeriod);
G_X = zeros(1,forecastPeriod);
G_A = zeros(1,forecastPeriod);
G_k = zeros(1,forecastPeriod);
G_T = zeros(1,forecastPeriod);

for j = 1:forecastPeriod
   
    G_C(1,j) = sum(G_C_temp(1,1:j));
    G_w(1,j) = sum(G_w_temp(1,1:j));
    G_X(1,j) = sum(G_X_temp(1,1:j));
    G_A(1,j) = sum(G_A_temp(1,1:j));
    G_k(1,j) = sum(G_k_temp(1,1:j));
    G_T(1,j) = sum(G_T_temp(1,1:j));
    
end

clear j

% now initialize learning parameters

phiR  = zeros(3,timePeriod); % it is 3 because we have 3 variables (2 state variables + intercept)
phiW  = zeros(3,timePeriod);
phiPi = zeros(3,timePeriod);
phiT  = zeros(3,timePeriod);
phiK  = zeros(3,timePeriod);

% now initialize coefs as small random numbers

phiR(:,1)  = 0.001*randn(3,1);
phiW(:,1)  = 0.001*randn(3,1);
phiPi(:,1) = 0.001*randn(3,1);
phiT(:,1)  = 0.001*randn(3,1);
phiK(:,1)  = 0.001*randn(3,1);
D_Old      = eye(3); % initialize D as an identity matrix (meaning there are no cors between variables

% compute Lst coefficient

LsCoef = (1/( ((eta-alpha)*Lst^(eta-alpha)) - ((1-alpha)/(alpha*cst))*(1/beta)*alphaY*(1-alpha)) )*((1-alpha)/(alpha*cst));

% do learning

for t = 2:timePeriod
    
    if t == 2
       
       % initialize at steady state
        
       allVariables(1,t-1)  = Lst;
       allVariables(2,t-1)  = kst;
       allVariables(3,t-1)  = Yst;
       allVariables(4,t-1)  = Invst;
       allVariables(5,t-1)  = cst;
       allVariables(6,t-1)  = rKst;
       allVariables(7,t-1)  = wst;
       allVariables(8,t-1)  = Rst;
       allVariables(9,t-1)  = profit_st;
       allVariables(10,t-1) = Tst;
       allVariables(11,t-1) = 1; % this is inflation
       allVariables(12,t-1) = Xst;
       
       % LEARNING ALGO, BEGIN
       
       zMat  = [ 1 allVariables(2,t-1) A_hat(1,t-1) ]'; % this is capital and productivity shock at T-1
       D_New = D_Old + gainParam * (zMat * zMat' - D_Old);     % update D matrix
       
       phiR(:,t)  = phiR(:,t-1)  + gainParam * D_New \ zMat * (allVariables(8,t)  - phiR(:,t-1)'  * zMat); % update learning 
       phiW(:,t)  = phiW(:,t-1)  + gainParam * D_New \ zMat * (allVariables(7,t)  - phiW(:,t-1)'  * zMat);
       phiPi(:,t) = phiPi(:,t-1) + gainParam * D_New \ zMat * (allVariables(11,t) - phiPi(:,t-1)' * zMat);
       phiT(:,t)  = phiT(:,t-1)  + gainParam * D_New \ zMat * (allVariables(10,t) - phiT(:,t-1)'  * zMat);
       phiK(:,t)  = phiK(:,t-1)  + gainParam * D_New \ zMat * (allVariables(2,t)  - phiK(:,t-1)'  * zMat);
       
       % LEARNING ALOG, END
       
       % calculate predictions based on estimates for one period ahead
       
       CGpredictedR(1,t-1)  = phiR(:,t)'  * zMat; % R
       CGpredictedW(1,t-1)  = phiW(:,t)'  * zMat; % Wage
       CGpredictedPi(1,t-1) = phiPi(:,t)' * zMat; % Inflation
       CGpredictedT(1,t-1)  = phiT(:,t)'  * zMat; % Transfers
       CGpredictedK(1,t-1)  = phiK(:,t)'  * zMat; % Capital
       
       % calculate predictions j steps ahead
       
       % initialize predictions
       
       RPrediction  = zeros(1,forecastPeriod);
       wPrediction  = zeros(1,forecastPeriod);
       PiPrediction = zeros(1,forecastPeriod);
       TPrediction  = zeros(1,forecastPeriod);
       KPrediction  = zeros(1,forecastPeriod);
    
       RPrediction(1,1)  = CGpredictedR(1,t-1);
       wPrediction(1,1)  = CGpredictedW(1,t-1);
       PiPrediction(1,1) = CGpredictedPi(1,t-1);
       TPrediction(1,1)  = CGpredictedT(1,t-1);
       KPrediction(1,1)  = CGpredictedK(1,t-1);
       
       % create zMatrix for computing predictions
       
       zMatPred = [ 1 CGpredictedK(1,t-1) A_hat(1,t-1) ]';
       
       
       for i = 2:forecastPeriod
           
            RPrediction(1,i)  = phiR(:,t)'  * zMatPred;
            wPrediction(1,i)  = phiW(:,t)'  * zMatPred;
            PiPrediction(1,i) = phiPi(:,t)' * zMatPred;
            TPrediction(1,i)  = phiT(:,t)'  * zMatPred;
            KPrediction(1,i)  = phiK(:,t)'  * zMatPred;
            
            zMatPred = [ 1 RPrediction(1,i-1) A_hat(1,t+i-1) ]';
            
       end
       
       
       % now compute big sums with those i steps ahead forecasts from above
       
       R_Pi_Sum   = zeros(1,forecastPeriod);
       R_R_Sum    = zeros(1,forecastPeriod);
       
       G_C_R_Sum  = zeros(1,forecastPeriod);
       G_C_Pi_Sum = zeros(1,forecastPeriod);
       G_W_W_Sum  = zeros(1,forecastPeriod);
       G_A_A_Sum  = zeros(1,forecastPeriod);
       G_K_Sum    = zeros(1,forecastPeriod);
       G_T_Sum    = zeros(1,forecastPeriod);
       Big_Sum_Pi = zeros(1,forecastPeriod);
       Big_Sum_R  = zeros(1,forecastPeriod);
       Sum_2_Pi   = zeros(1,forecastPeriod);
       
       for i = 1:forecastPeriod
          % is this correct? Let's see...
          
          G_C_R_Sum(1,i)  = G_C(1,forecastPeriod) * sum(RPrediction(1,1:i));
          G_C_Pi_Sum(1,i) = G_C(1,forecastPeriod) * sum(PiPrediction(1,1:i));
          G_W_W_Sum(1,i)  = G_w(1,i) * sum(wPrediction(1,1:i));
          G_A_A_Sum(1,i)  = G_A(1,i) * A_hat(1,t+i);
          G_K_Sum(1,i)    = G_k(1,i) * sum(KPrediction(1,1:i));
          G_T_Sum(1,i)    = G_T(1,i) * sum(TPrediction(1,1:i));
          
          Big_Sum_Pi(1,i) = BigLambda *( ((1/Rst)^i)*sum(PiPrediction(1,1:i)) );
          Big_Sum_R(1,i)  = BigLambda *( ((1/Rst)^i)*sum(RPrediction(1,1:i)) );
          
          Sum_2_Pi(1,i)   = G_X(1,i) * (PiPrediction(1,t) - allVariables(9,t-1));
          
       end
       
       % capital
       
       allVariables(2,t) = (1-delta)*allVariables(2,t-1) + allVariables(4,t-1);
       
       % inflation, Phillips curve
      
       allVariables(11,t) = (1/beta)*(allVariables(11,t-1)+(1/theta)*((1-theta*beta)*(1-theta))*allVariables(12,t-1));
       
       % consumption
       allVariables(5,t) = (-1/(C_0+G_C(1,forecastPeriod))) * (C_w*CGpredictedW(1,t-1) + (theta/((1-theta*beta)*(1-theta)))*C_X*(beta*CGpredictedPi(1,t-1)-allVariables(11,t)) + ...
                           C_A*A_hat(1,t) + C_k *allVariables(2,t) + Tst*CGpredictedT(1,t-1) + Rst*kst*(CGpredictedR(1,t-1)+allVariables(2,t-1)-allVariables(9,t-1)) + ...
                           sum(Big_Sum_Pi) - sum(Big_Sum_R) + sum(G_C_R_Sum) - sum(G_C_Pi_Sum) + sum(G_W_W_Sum) + (theta/((1-theta*beta)*(1-theta)))*sum(Sum_2_Pi) +...
                           sum(G_A_A_Sum) + sum(G_K_Sum) + sum(G_T_Sum));
      
      
       
      % labour
      
      allVariables(1,t) = LsCoef*((1/beta)*(alphaPi*allVariables(11,t)+alphaY*(A_hat(1,t) + alpha*CGpredictedK(1,t-1))) - allVariables(11,t) - ((1/beta) - 1 - delta)*allVariables(5,t));
      
      % output
      
      allVariables(3,t) = A_hat(1,t) + alpha*(allVariables(2,t)) + (1-alpha)*(allVariables(1,t));
      
      % wage
      
      allVariables(7,t) = allVariables(5,t) + (1-eta)*allVariables(1,t);
      
      % markup Xst
      
      allVariables(12,t) = A_hat(1,t) + alpha*allVariables(2,t) - alpha*allVariables(1,t) - allVariables(7,t);
      
      % solve for rK
      
      allVariables(6,t) = (alpha-1)*allVariables(2,t) - allVariables(12,t) + (1-alpha)*allVariables(1,t);
      
      % solve Investment
      
      allVariables(4,t) = allVariables(3,t)*(Yst/Invst) - (cst/Invst)*allVariables(5,t);
      
      % profit
      
      allVariables(9,t) = (Yst*(1-(1/Xst)) *allVariables(3,t) + (Yst/Xst)*allVariables(12,t))/profit_st;
      
      % R interest rate, Taylor rule
      
      allVariables(8,t) = alphaPi*allVariables(11,t) + alphaY*allVariables(3,t);
      
      % transfers
      
      allVariables(10,t) = (cst*allVariables(5,t) + kst*allVariables(2,t) - wst*Lst*allVariables(7,t) - wst*Lst*allVariables(1,t)-Rst*kst*allVariables(8,t) - Rst*kst*allVariables(2,t) + ...
                            Rst*kst*allVariables(11,t) - profit_st*allVariables(8,t))/Tst;
                        
    else
        
        
       zMat  = [ 1 allVariables(2,t-1) A_hat(1,t-1) ]'; % this is capital and productivity shock at T-1
       D_New = D_Old + gainParam * (zMat * zMat' - D_Old);     % update D matrix
        
       phiR(:,t)  = phiR(:,t-1)  + gainParam * D_New \ zMat * (allVariables(8,t)  - phiR(:,t-1)'  * zMat); % update learning
       phiW(:,t)  = phiW(:,t-1)  + gainParam * D_New \ zMat * (allVariables(7,t)  - phiW(:,t-1)'  * zMat);
       phiPi(:,t) = phiPi(:,t-1) + gainParam * D_New \ zMat * (allVariables(11,t) - phiPi(:,t-1)' * zMat);
       phiT(:,t)  = phiT(:,t-1)  + gainParam * D_New \ zMat * (allVariables(10,t) - phiT(:,t-1)'  * zMat);
       phiK(:,t)  = phiK(:,t-1)  + gainParam * D_New \ zMat * (allVariables(2,t)  - phiK(:,t-1)'  * zMat);
       
       CGpredictedR(1,t-1)  = phiR(:,t)'  * zMat; % R
       CGpredictedW(1,t-1)  = phiW(:,t)'  * zMat; % Wage
       CGpredictedPi(1,t-1) = phiPi(:,t)' * zMat; % Inflation
       CGpredictedT(1,t-1)  = phiT(:,t)'  * zMat; % Transfers
       CGpredictedK(1,t-1)  = phiK(:,t)'  * zMat; % Capital
       
       RPrediction  = zeros(1,forecastPeriod);
       wPrediction  = zeros(1,forecastPeriod);
       PiPrediction = zeros(1,forecastPeriod);
       TPrediction  = zeros(1,forecastPeriod);
       KPrediction  = zeros(1,forecastPeriod);
       
       RPrediction(1,1)  = CGpredictedR(1,t-1);
       wPrediction(1,1)  = CGpredictedW(1,t-1);
       PiPrediction(1,1) = CGpredictedPi(1,t-1);
       TPrediction(1,1)  = CGpredictedT(1,t-1);
       KPrediction(1,1)  = CGpredictedK(1,t-1);
       
       zMatPred = [ 1 CGpredictedK(1,t-1) A_hat(1,t-1) ]';
       
       clear i
       
       for i = 2:forecastPeriod
           
            RPrediction(1,i)  = phiR(:,t)'  * zMatPred;
            wPrediction(1,i)  = phiW(:,t)'  * zMatPred;
            PiPrediction(1,i) = phiPi(:,t)' * zMatPred;
            TPrediction(1,i)  = phiT(:,t)'  * zMatPred;
            KPrediction(1,i)  = phiK(:,t)'  * zMatPred;
            
            zMatPred = [ 1 RPrediction(1,i-1) A_hat(1,t+i-1) ]';
            
       end
       
       R_Pi_Sum   = zeros(1,forecastPeriod);
       R_R_Sum    = zeros(1,forecastPeriod);
       
       G_C_R_Sum  = zeros(1,forecastPeriod);
       G_C_Pi_Sum = zeros(1,forecastPeriod);
       G_W_W_Sum  = zeros(1,forecastPeriod);
       G_A_A_Sum  = zeros(1,forecastPeriod);
       G_K_Sum    = zeros(1,forecastPeriod);
       G_T_Sum    = zeros(1,forecastPeriod);
       Big_Sum_Pi = zeros(1,forecastPeriod);
       Big_Sum_R  = zeros(1,forecastPeriod);
       Sum_2_Pi   = zeros(1,forecastPeriod);
       
       clear i
       
       for i = 1:forecastPeriod
           
          % is this correct? Let's see...
          
          G_C_R_Sum(1,i)  = G_C(1,forecastPeriod) * sum(RPrediction(1,1:i));
          G_C_Pi_Sum(1,i) = G_C(1,forecastPeriod) * sum(PiPrediction(1,1:i));
          G_W_W_Sum(1,i)  = G_w(1,i) * sum(wPrediction(1,1:i));
          G_A_A_Sum(1,i)  = G_A(1,i) * A_hat(1,t+i);
          G_K_Sum(1,i)    = G_k(1,i) * sum(KPrediction(1,1:i));
          G_T_Sum(1,i)    = G_T(1,i) * sum(TPrediction(1,1:i));
          
          Big_Sum_Pi(1,i) = BigLambda *( ((1/Rst)^i)*sum(PiPrediction(1,1:i)) );
          Big_Sum_R(1,i)  = BigLambda *( ((1/Rst)^i)*sum(RPrediction(1,1:i)) );
          
          Sum_2_Pi(1,i)   = G_X(1,i) * (PiPrediction(1,t) - allVariables(9,t-1));
          
       end
       
        % capital
       
       allVariables(2,t) = (1-delta)*allVariables(2,t-1) + allVariables(4,t-1);
       
       % inflation, Phillips curve
      
       allVariables(11,t) = (1/beta)*(allVariables(11,t-1)+(1/theta)*((1-theta*beta)*(1-theta))*allVariables(12,t-1));
       
       % consumption
       allVariables(5,t) = (-1/(C_0+G_C(1,forecastPeriod))) * (C_w*CGpredictedW(1,t-1) + (theta/((1-theta*beta)*(1-theta)))*C_X*(beta*CGpredictedPi(1,t-1)-allVariables(11,t)) + ...
                           C_A*A_hat(1,t) + C_k *allVariables(2,t) + Tst*CGpredictedT(1,t-1) + Rst*kst*(CGpredictedR(1,t-1)+allVariables(2,t-1)-allVariables(9,t-1)) + ...
                           sum(Big_Sum_Pi) - sum(Big_Sum_R) + sum(G_C_R_Sum) - sum(G_C_Pi_Sum) + sum(G_W_W_Sum) + (theta/((1-theta*beta)*(1-theta)))*sum(Sum_2_Pi) +...
                           sum(G_A_A_Sum) + sum(G_K_Sum) + sum(G_T_Sum));
      
      
       
      % labour
      
      allVariables(1,t) = LsCoef*((1/beta)*(alphaPi*allVariables(11,t)+alphaY*(A_hat(1,t) + alpha*CGpredictedK(1,t-1))) - allVariables(11,t) - ((1/beta) - 1 - delta)*allVariables(5,t));
      
      % output
      
      allVariables(3,t) = A_hat(1,t) + alpha*(allVariables(2,t)) + (1-alpha)*(allVariables(1,t));
      
      % wage
      
      allVariables(7,t) = allVariables(5,t) + (1-eta)*allVariables(1,t);
      
      % markup Xst
      
      allVariables(12,t) = A_hat(1,t) + alpha*allVariables(2,t) - alpha*allVariables(1,t) - allVariables(7,t);
      
      % solve for rK
      
      allVariables(6,t) = (alpha-1)*allVariables(2,t) - allVariables(12,t) + (1-alpha)*allVariables(1,t);
      
      % solve Investment
      
      allVariables(4,t) = allVariables(3,t)*(Yst/Invst) - (cst/Invst)*allVariables(5,t);
      
      % profit
      
      allVariables(9,t) = (Yst*(1-(1/Xst)) *allVariables(3,t) + (Yst/Xst)*allVariables(12,t))/profit_st;
      
      % R interest rate, Taylor rule
      
      allVariables(8,t) = alphaPi*allVariables(11,t) + alphaY*allVariables(3,t);
      
      % transfers
      
      allVariables(10,t) = (cst*allVariables(5,t) + kst*allVariables(2,t) - wst*Lst*allVariables(7,t) - wst*Lst*allVariables(1,t)-Rst*kst*allVariables(8,t) - Rst*kst*allVariables(2,t) + ...
                            Rst*kst*allVariables(11,t) - profit_st*allVariables(8,t))/Tst;
       
       
    end
    
    
end







































