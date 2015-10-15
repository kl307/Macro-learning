// This is dynare file for New Keynesian model with learning. The main reference for this code is Gali's book (2nd chapter)

var R w c L rK k Y Inv X infl A;

predetermined_variables k;

varexo shock;

parameters eta beta epsilon alpha delta rho theta alphaPI alphaY Lst kst Yst Invst cst rKst wst Rst;

alphaPI = 1.3;
alphaY  = 0.5;

// variables are:
// R    - nominal interest rate
// w    - wage paid to household
// c    - household consumption
// L    - labour
// rK   - real interest (on capital)
// k    - capital
// Y    - output
// Inv  - investment
// X    - price markup (sticky prices)
// infl - final good inflation
// A    - autoregressive productivity shock

// parameters are:
// eta     - labor supply elasticty
// beta    - discount factor
// epsilon - price markup
// alpha   - production function parameter
// delta   - capital depreciation
// rho     - autoregressive parameter in the shock
// alphaPI - taylor rule parameter
// alphaY  - taylor rule parameter

model(linear); 

// 1 eqn. Phillips Curve

beta*infl(+1) = infl + ((1-beta*theta)*(1-theta)/theta) * X;

// 2 eqn. Shock Equation

A = rho*A(-1) + shock;

// 3 eqn. Household Euler equation 

-c = -c(+1) + R(+1) - infl(+1);

// 4 eqn. Labour supply

w = c + (eta-1)*L;

// 5 eqn. Real Interest Rate

rK = (alpha-1)*k - X + (1-alpha)*L + A;

// 6 eqn. Some Identity

R*Rst = rK*rKst + Rst*infl;

// 7 eqn. Labour Demand

w = A + alpha*k - alpha*L - X;

// 8 eqn. Production Function

Y = alpha*k + (1-alpha)*L;

// 9 eqn. Capital Accumulation

k(+1) = (1-delta)*k + delta*Inv;

// 10 eqn. Production and Investment

Y = (cst/Yst)*c + (Invst/Yst)*Inv;

// 11 eqn. Taylor Rule

R = alphaPI*infl + alphaY*Y;

end;

shocks;

var shock;

stderr 0.01;

end;

resid(1);

steady;

check;

stoch_simul(irf=10,order=1);

conditional_variance_decomposition=1;

















