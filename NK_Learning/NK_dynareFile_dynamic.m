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

residual = zeros(11, 1);
lhs =params(2)*y(15);
rhs =y(12)+(1-params(2)*params(7))*(1-params(7))/params(7)*y(11);
residual(1)= lhs-rhs;
lhs =y(13);
rhs =params(6)*y(2)+x(it_, 1);
residual(2)= lhs-rhs;
lhs =y(5);
rhs =y(1)+y(3)-y(12);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(5)+(params(1)-1)*y(6);
residual(4)= lhs-rhs;
lhs =y(7);
rhs =(params(4)-1)*y(8)-y(11)+y(6)*(1-params(4));
residual(5)= lhs-rhs;
lhs =y(3)*params(17);
rhs =y(7)*params(15)-y(12)*params(17);
residual(6)= lhs-rhs;
lhs =y(4);
rhs =y(13)+params(4)*y(8)+y(6)*params(4)-y(11);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(6)*(1-params(4))+params(4)*y(8);
residual(8)= lhs-rhs;
lhs =y(14);
rhs =y(8)*(1-params(5))+params(5)*y(10);
residual(9)= lhs-rhs;
lhs =y(9);
rhs =y(5)*params(14)/params(12)+y(10)*params(13)/params(12);
residual(10)= lhs-rhs;
lhs =y(3);
rhs =y(12)*params(8)+y(9)*params(9);
residual(11)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(11, 16);

  %
  % Jacobian matrix
  %

  g1(1,11)=(-((1-params(2)*params(7))*(1-params(7))/params(7)));
  g1(1,12)=(-1);
  g1(1,15)=params(2);
  g1(2,2)=(-params(6));
  g1(2,13)=1;
  g1(2,16)=(-1);
  g1(3,3)=(-1);
  g1(3,1)=(-1);
  g1(3,5)=1;
  g1(3,12)=1;
  g1(4,4)=1;
  g1(4,5)=(-1);
  g1(4,6)=(-(params(1)-1));
  g1(5,6)=(-(1-params(4)));
  g1(5,7)=1;
  g1(5,8)=(-(params(4)-1));
  g1(5,11)=1;
  g1(6,3)=params(17);
  g1(6,7)=(-params(15));
  g1(6,12)=params(17);
  g1(7,4)=1;
  g1(7,6)=(-params(4));
  g1(7,8)=(-params(4));
  g1(7,11)=1;
  g1(7,13)=(-1);
  g1(8,6)=(-(1-params(4)));
  g1(8,8)=(-params(4));
  g1(8,9)=1;
  g1(9,8)=(-(1-params(5)));
  g1(9,14)=1;
  g1(9,10)=(-params(5));
  g1(10,5)=(-(params(14)/params(12)));
  g1(10,9)=1;
  g1(10,10)=(-(params(13)/params(12)));
  g1(11,3)=1;
  g1(11,9)=(-params(9));
  g1(11,12)=(-params(8));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,256);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],11,4096);
end
end
