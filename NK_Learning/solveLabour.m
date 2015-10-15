function LabourSteadyState = solveLabour(Lst)

    % Purpuse of this function:
    % evaluate expression at zero to obtain steady state value of labour
    % we load struct named 'paramValues' which has all the parameter values
    % needed to evaluate the expression
    
    load paramValues
    
    eta     = paramValues.eta;
    beta    = paramValues.beta;
    alpha   = paramValues.alpha;
    delta   = paramValues.delta;
    Xst     = paramValues.Xst;
    
    clear paramValues
    
    LabourSteadyState = ((((1-alpha)/Xst)*((Xst/alpha)*((1/beta)-1+delta))^(alpha/(alpha-1)))*Lst) - ...
       (Lst^(eta-1))*(Lst*((Xst/alpha) *((1/beta)-1+delta))^(alpha/(alpha-1)) - delta*Lst*((Xst/alpha)-1+delta)^(1/(alpha-1)));
   

end