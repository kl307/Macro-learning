function output = find_h(x)

% This function solves the steady state value of the housing ratio
% distributed among households and entrepreneurs.

global epsilon v gamma beta paramj m_star_st

%find the steady state value of h

output= (1-gamma-m_star_st*beta*(1-gamma/beta))*paramj*(1+(1-v)*(epsilon-1))*x/(1-beta)/((epsilon-1)-(epsilon-1)*m_star_st*paramj*x)-gamma*v;

end

