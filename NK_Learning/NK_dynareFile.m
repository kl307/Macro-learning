%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'NK_dynareFile';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('NK_dynareFile.log');
M_.exo_names = 'shock';
M_.exo_names_tex = 'shock';
M_.exo_names_long = 'shock';
M_.endo_names = 'R';
M_.endo_names_tex = 'R';
M_.endo_names_long = 'R';
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'rK');
M_.endo_names_tex = char(M_.endo_names_tex, 'rK');
M_.endo_names_long = char(M_.endo_names_long, 'rK');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'Y');
M_.endo_names = char(M_.endo_names, 'Inv');
M_.endo_names_tex = char(M_.endo_names_tex, 'Inv');
M_.endo_names_long = char(M_.endo_names_long, 'Inv');
M_.endo_names = char(M_.endo_names, 'X');
M_.endo_names_tex = char(M_.endo_names_tex, 'X');
M_.endo_names_long = char(M_.endo_names_long, 'X');
M_.endo_names = char(M_.endo_names, 'infl');
M_.endo_names_tex = char(M_.endo_names_tex, 'infl');
M_.endo_names_long = char(M_.endo_names_long, 'infl');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
M_.endo_names = char(M_.endo_names, 'T');
M_.endo_names_tex = char(M_.endo_names_tex, 'T');
M_.endo_names_long = char(M_.endo_names_long, 'T');
M_.endo_names = char(M_.endo_names, 'profit');
M_.endo_names_tex = char(M_.endo_names_tex, 'profit');
M_.endo_names_long = char(M_.endo_names_long, 'profit');
M_.param_names = 'eta';
M_.param_names_tex = 'eta';
M_.param_names_long = 'eta';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'alphaPI');
M_.param_names_tex = char(M_.param_names_tex, 'alphaPI');
M_.param_names_long = char(M_.param_names_long, 'alphaPI');
M_.param_names = char(M_.param_names, 'alphaY');
M_.param_names_tex = char(M_.param_names_tex, 'alphaY');
M_.param_names_long = char(M_.param_names_long, 'alphaY');
M_.param_names = char(M_.param_names, 'Lst');
M_.param_names_tex = char(M_.param_names_tex, 'Lst');
M_.param_names_long = char(M_.param_names_long, 'Lst');
M_.param_names = char(M_.param_names, 'kst');
M_.param_names_tex = char(M_.param_names_tex, 'kst');
M_.param_names_long = char(M_.param_names_long, 'kst');
M_.param_names = char(M_.param_names, 'Yst');
M_.param_names_tex = char(M_.param_names_tex, 'Yst');
M_.param_names_long = char(M_.param_names_long, 'Yst');
M_.param_names = char(M_.param_names, 'Invst');
M_.param_names_tex = char(M_.param_names_tex, 'Invst');
M_.param_names_long = char(M_.param_names_long, 'Invst');
M_.param_names = char(M_.param_names, 'cst');
M_.param_names_tex = char(M_.param_names_tex, 'cst');
M_.param_names_long = char(M_.param_names_long, 'cst');
M_.param_names = char(M_.param_names, 'rKst');
M_.param_names_tex = char(M_.param_names_tex, 'rKst');
M_.param_names_long = char(M_.param_names_long, 'rKst');
M_.param_names = char(M_.param_names, 'wst');
M_.param_names_tex = char(M_.param_names_tex, 'wst');
M_.param_names_long = char(M_.param_names_long, 'wst');
M_.param_names = char(M_.param_names, 'Rst');
M_.param_names_tex = char(M_.param_names_tex, 'Rst');
M_.param_names_long = char(M_.param_names_long, 'Rst');
M_.param_names = char(M_.param_names, 'profit_st');
M_.param_names_tex = char(M_.param_names_tex, 'profit\_st');
M_.param_names_long = char(M_.param_names_long, 'profit_st');
M_.param_names = char(M_.param_names, 'Xst');
M_.param_names_tex = char(M_.param_names_tex, 'Xst');
M_.param_names_long = char(M_.param_names_long, 'Xst');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 13;
M_.param_nbr = 19;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('NK_dynareFile_static');
erase_compiled_function('NK_dynareFile_dynamic');
M_.lead_lag_incidence = [
 0 3 0;
 0 4 0;
 0 5 16;
 0 6 0;
 0 7 0;
 1 8 0;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 17;
 2 13 0;
 0 14 0;
 0 15 0;]';
M_.nstatic = 9;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(19, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 50;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 8 ) = 1.3;
alphaPI = M_.params( 8 );
M_.params( 9 ) = 0.5;
alphaY = M_.params( 9 );
M_.params( 19 ) = M_.params(3)/(M_.params(3)-1);
Xst = M_.params( 19 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.01)^2;
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 40;
options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);
conditional_variance_decomposition=1;
save('NK_dynareFile_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('NK_dynareFile_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('NK_dynareFile_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('NK_dynareFile_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('NK_dynareFile_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
