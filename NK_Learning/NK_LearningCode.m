% this is learning code for NK Model

% 1. Run this file to obtain steady state solution and parameter values

steadyStateSolution;

clear paramValues Linit

timePeriod      = 40;   % time period of model
forecastPeriod  = 100;  % forecasting period
gainParam       = 0.04; % constant gain parameter
numVars         = 12;   % number of variables in the model
allVariablesT   = zeros(numVars,timePeriod); % all variables at time T
allVariablesLag = zeros(numVars,timePeriod); % all variables at time T-1

% create vector of shocks

A_hat      = zeros(1,timePeriod+1);
A_hat(1,1) = 1; % initialize A_hat at 1

for j = 2:timePeriod+1
   
    A_hat(1,j) = 0.001*randn;
    
end

% initialize the productivity shock

epsA = 0.01*randn;

allVariablesT(12,:)   = A_hat(2:end);
allVariablesLag(12,:) = A_hat(1:end-1);


% define some parameters for consumption function

C_w = wst*( (eta/(eta-1))*(wst^(1/(eta-1)))*(cst^(-1/(eta-1)))+((1-alpha)/(eta-1))*(1-(1/Xst))*(kst^alpha)*(cst^(-(1-alpha)/(eta-1)))*(wst^((1-alpha)/(eta-1)-1)) );
C_0 = cst*( -(eta/(eta-1))*(wst^(eta/(eta-1)))*(cst^(-eta/(eta-1)))-((1-alpha)/(eta-1))*(1-(1/Xst))*(kst^alpha)*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1)-1)) -1 );
C_X = (1/Xst)*( (kst^alpha)*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1))) );
C_A = (1-(1/Xst))*(kst^alpha)*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1)));
C_k = alpha*(1-(1/Xst))*(kst^(1-alpha))*(wst^((1-alpha)/(eta-1)))*(cst^(-(1-alpha)/(eta-1)));

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

% now initialize learning parameters

phiR  = zeros(3,timePeriod); % it is 3 because we have 3 variables (2 state variables + intercept)
phiW  = zeros(3,timePeriod);
phiPi = zeros(3,timePeriod);
phiT  = zeros(3,timePeriod);

% now initialize coefs as small random numbers

phiR(:,1)  = 0.001*randn(3,1);
phiW(:,1)  = 0.001*randn(3,1);
phiPi(:,1) = 0.001*randn(3,1);
phiT(:,1)  = 0.001*randn(3,1);
D_old      = eye(3); % initialize D as an identity matrix (meaning there are no cors between variables

% do learning

for t = 2:timePeriod
    
    if t == 2
       
       % initialize at steady state
        
       allVariablesLag(1,t-1)  = Lst;
       allVariablesLag(2,t-1)  = kst;
       allVariablesLag(3,t-1)  = Yst;
       allVariablesLag(4,t-1)  = Invst;
       allVariablesLag(5,t-1)  = cst;
       allVariablesLag(6,t-1)  = rKst;
       allVariablesLag(7,t-1)  = wst;
       allVariablesLag(8,t-1)  = Rst;
       allVariablesLag(9,t-1)  = profit_st;
       allVariablesLag(10,t-1) = Tst;
       allVariablesLag(11,t-1) = 1; % this is inflation
       
       - = ( 1/(C_0 + sum(G_C)) )
       
       
    else
        
        %here we write how economy evolves
        
       allVariablesT(1,t)  = 
       allVariablesT(3,t)  = 
       allVariablesT(4,t)  = 
       allVariablesT(5,t)  = 
       allVariablesT(6,t)  = 
       allVariablesT(7,t)  = 
       allVariablesT(8,t)  = 
       allVariablesT(9,t)  = 
       allVariablesT(10,t) = 
       
        
    end
    
    zMat  = [ 1 allVariablesLag(2,t) allVariablesLag(11,t) ]; % this is capital and productivity shock at T-1
    D_New = D_Old + gainParam * (zMat * zMat' - D_Old);    % update D matrix
    
    % update parameters of PLM using constant gain learning
    
    phiR(:,t)  = phiR(:,t-1)  + gainParam * D_New \ zMat * (allVariablesT(8,t)  - phiR(:,t-1)'  * zMat);
    phiW(:,t)  = phiW(:,t-1)  + gainParam * D_New \ zMat * (allVariablesT(7,t)  - phiW(:,t-1)'  * zMat);
    phiPi(:,t) = phiPi(:,t-1) + gainParam * D_New \ zMat * (allVariablesT(11,t) - phiPi(:,t-1)' * zMat);
    phiT(:,t)  = phiT(:,t-1)  + gainParam * D_New \ zMat * (allVariables(10,t)  - phiT(:,t-1)'  * zMat);
        
    % update D and Z and calculate the predicted (expected) variables
    
    D_Old = D_New;
    
    zMat  = [ 1 allVariablesT(2,t) allVariablesT(11,t) ]; % this is capital and productivity shock at T-1
    
    % now we predict all the variables using the coefficients we obtained
    
    predictedR  = phiR(:,t)'  * zMat;
    predictedW  = phiW(:,t)'  * zMat;
    predictedPi = phiPi(:,t)' * zMat;
    predictedT  = phiT(:,t)'  * zMat;
    
    % now we have to compute forecasts some steps ahead using the same
    % coefficients we obtained
    
    % initialize first
    
    RPrediction  = zeros(1,forecastPeriod);
    wPrediction  = zeros(1,forecastPeriod);
    PiPrediction = zeros(1,forecastPeriod);
    TPrediction  = zeros(1,forecastPeriod);
    
    % initialize predictions
    
    RPrediction(1,1)  = predictedR;
    wPrediction(1,1)  = predictedW;
    PiPrediction(1,1) = predictedPi;
    TPrediction(1,1)  = predictedT;
    
    for i = 2:forecastPeriod
        
        RPrediction(1,i)  = phiR(:,t)'  * RPrediction(1,i-1);
        wPrediction(1,i)  = phiW(:,t)'  * wPrediction(1,i-1);
        PiPrediction(1,i) = phiPi(:,t)' * PiPrediction(1,i-1);
        TPrediction(1,i)  = phiT(:,t)'  * TPrediction(1,i-1);
        
    end
    
    % now we can obtain consumption knowing all the other variables
    
    
    
    
    
    
    
    
end







































