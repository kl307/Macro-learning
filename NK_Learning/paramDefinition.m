% this file defines parameter values

eta     = 2;
beta    = 0.98;
epsilon = 6;
alpha   = 0.3;
delta   = 0.025;
theta   = 0.25;
rho     = 0.9;
alphaPi = 0.5;
alphaY  = 0.45;
    
paramValues = struct('eta',eta,'beta',beta,'epsilon',epsilon,'alpha',alpha,'delta',delta,'Xst',epsilon/(epsilon--1),'theta',theta,'rho',rho,...
                     'alphaPI',alphaPi,'alphaY',alphaY);
    
clear eta beta epsilon alpha delta Xst theta rho
save paramValues paramValues
