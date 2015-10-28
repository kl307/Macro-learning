function [residual, g1, g2, g3] = NK_dynareFile_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(14, 1);
lhs =params(2)*y(20);
rhs =y(13)+(1-params(2)*params(7))*(1-params(7))/params(7)*y(12);
residual(1)= lhs-rhs;
lhs =y(14);
rhs =params(6)*y(2)+x(it_, 1);
residual(2)= lhs-rhs;
lhs =(-y(6));
rhs =(-y(19))+y(18)-y(20);
residual(3)= lhs-rhs;
lhs =y(5);
rhs =y(6)+(params(1)-1)*y(7);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =y(14)+y(7)*(1-params(4))+(params(4)-1)*y(1)-y(12);
residual(5)= lhs-rhs;
lhs =y(4)*params(17);
rhs =y(8)*params(15)+y(13)*params(17);
residual(6)= lhs-rhs;
lhs =y(5);
rhs =y(14)+params(4)*y(1)-y(7)*params(4)-y(12);
residual(7)= lhs-rhs;
lhs =y(10);
rhs =y(7)*(1-params(4))+params(4)*y(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =params(5)*y(11)+(1-params(5))*y(1);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(6)*params(14)/params(12)+y(11)*params(13)/params(12);
residual(10)= lhs-rhs;
lhs =y(4);
rhs =y(13)*params(8)+y(10)*params(9);
residual(11)= lhs-rhs;
lhs =y(16)*params(18);
rhs =y(10)*params(12)*(1-1/params(19))+y(12)*params(12)/params(19);
residual(12)= lhs-rhs;
lhs =y(15);
rhs =y(13)*params(17)*params(11)+y(6)*params(14)+params(11)*y(1)-y(7)*params(16)*params(10)-y(5)*params(16)*params(10)-y(16)*params(18)-y(4)*params(17)*params(11)-params(17)*params(11)*y(3);
residual(13)= lhs-rhs;
lhs =y(17);
rhs =y(1);
residual(14)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(14, 21);

  %
  % Jacobian matrix
  %

  g1(1,12)=(-((1-params(2)*params(7))*(1-params(7))/params(7)));
  g1(1,13)=(-1);
  g1(1,20)=params(2);
  g1(2,2)=(-params(6));
  g1(2,14)=1;
  g1(2,21)=(-1);
  g1(3,18)=(-1);
  g1(3,6)=(-1);
  g1(3,19)=1;
  g1(3,20)=1;
  g1(4,5)=1;
  g1(4,6)=(-1);
  g1(4,7)=(-(params(1)-1));
  g1(5,7)=(-(1-params(4)));
  g1(5,8)=1;
  g1(5,1)=(-(params(4)-1));
  g1(5,12)=1;
  g1(5,14)=(-1);
  g1(6,4)=params(17);
  g1(6,8)=(-params(15));
  g1(6,13)=(-params(17));
  g1(7,5)=1;
  g1(7,7)=params(4);
  g1(7,1)=(-params(4));
  g1(7,12)=1;
  g1(7,14)=(-1);
  g1(8,7)=(-(1-params(4)));
  g1(8,1)=(-params(4));
  g1(8,10)=1;
  g1(9,1)=(-(1-params(5)));
  g1(9,9)=1;
  g1(9,11)=(-params(5));
  g1(10,6)=(-(params(14)/params(12)));
  g1(10,10)=1;
  g1(10,11)=(-(params(13)/params(12)));
  g1(11,4)=1;
  g1(11,10)=(-params(9));
  g1(11,13)=(-params(8));
  g1(12,10)=(-(params(12)*(1-1/params(19))));
  g1(12,12)=(-(params(12)/params(19)));
  g1(12,16)=params(18);
  g1(13,4)=params(17)*params(11);
  g1(13,5)=params(16)*params(10);
  g1(13,6)=(-params(14));
  g1(13,7)=params(16)*params(10);
  g1(13,1)=(-params(11));
  g1(13,13)=(-(params(17)*params(11)));
  g1(13,15)=1;
  g1(13,16)=params(18);
  g1(13,3)=params(17)*params(11);
  g1(14,1)=(-1);
  g1(14,17)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,441);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],14,9261);
end
end
