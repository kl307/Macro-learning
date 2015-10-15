function [residual, g1, g2] = NK_dynareFile_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 11, 1);

%
% Model equations
%

lhs =params(2)*y(10);
rhs =y(10)+(1-params(2)*params(7))*(1-params(7))/params(7)*y(9);
residual(1)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(6)+x(1);
residual(2)= lhs-rhs;
lhs =y(3);
rhs =y(3)+y(1)-y(10);
residual(3)= lhs-rhs;
lhs =y(2);
rhs =y(3)+(params(1)-1)*y(4);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =(params(4)-1)*y(6)-y(9)+y(4)*(1-params(4));
residual(5)= lhs-rhs;
lhs =y(1)*params(17);
rhs =y(5)*params(15)-y(10)*params(17);
residual(6)= lhs-rhs;
lhs =y(2);
rhs =y(11)+params(4)*y(6)+y(4)*params(4)-y(9);
residual(7)= lhs-rhs;
lhs =y(7);
rhs =y(4)*(1-params(4))+params(4)*y(6);
residual(8)= lhs-rhs;
lhs =y(6);
rhs =y(6)*(1-params(5))+params(5)*y(8);
residual(9)= lhs-rhs;
lhs =y(7);
rhs =y(3)*params(14)/params(12)+y(8)*params(13)/params(12);
residual(10)= lhs-rhs;
lhs =y(1);
rhs =y(10)*params(8)+y(7)*params(9);
residual(11)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(11, 11);

  %
  % Jacobian matrix
  %

  g1(1,9)=(-((1-params(2)*params(7))*(1-params(7))/params(7)));
  g1(1,10)=params(2)-1;
  g1(2,11)=1-params(6);
  g1(3,1)=(-1);
  g1(3,10)=1;
  g1(4,2)=1;
  g1(4,3)=(-1);
  g1(4,4)=(-(params(1)-1));
  g1(5,4)=(-(1-params(4)));
  g1(5,5)=1;
  g1(5,6)=(-(params(4)-1));
  g1(5,9)=1;
  g1(6,1)=params(17);
  g1(6,5)=(-params(15));
  g1(6,10)=params(17);
  g1(7,2)=1;
  g1(7,4)=(-params(4));
  g1(7,6)=(-params(4));
  g1(7,9)=1;
  g1(7,11)=(-1);
  g1(8,4)=(-(1-params(4)));
  g1(8,6)=(-params(4));
  g1(8,7)=1;
  g1(9,6)=1-(1-params(5));
  g1(9,8)=(-params(5));
  g1(10,3)=(-(params(14)/params(12)));
  g1(10,7)=1;
  g1(10,8)=(-(params(13)/params(12)));
  g1(11,1)=1;
  g1(11,7)=(-params(9));
  g1(11,10)=(-params(8));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,121);
end
end
