%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'nkm_3eq';
M_.dynare_version = '6.3';
oo_.dynare_version = '6.3';
options_.dynare_version = '6.3';
%
% Some global variables initialization
%
global_initialization;
options_.nograph = true;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'u_g'};
M_.exo_names_tex(1) = {'u\_g'};
M_.exo_names_long(1) = {'u_g'};
M_.endo_names = cell(5,1);
M_.endo_names_tex = cell(5,1);
M_.endo_names_long = cell(5,1);
M_.endo_names(1) = {'pi'};
M_.endo_names_tex(1) = {'pi'};
M_.endo_names_long(1) = {'pi'};
M_.endo_names(2) = {'ygap'};
M_.endo_names_tex(2) = {'ygap'};
M_.endo_names_long(2) = {'ygap'};
M_.endo_names(3) = {'i'};
M_.endo_names_tex(3) = {'i'};
M_.endo_names_long(3) = {'i'};
M_.endo_names(4) = {'rn'};
M_.endo_names_tex(4) = {'rn'};
M_.endo_names_long(4) = {'rn'};
M_.endo_names(5) = {'g'};
M_.endo_names_tex(5) = {'g'};
M_.endo_names_long(5) = {'g'};
M_.endo_partitions = struct();
M_.param_names = cell(10,1);
M_.param_names_tex = cell(10,1);
M_.param_names_long = cell(10,1);
M_.param_names(1) = {'gamma'};
M_.param_names_tex(1) = {'gamma'};
M_.param_names_long(1) = {'gamma'};
M_.param_names(2) = {'phi'};
M_.param_names_tex(2) = {'phi'};
M_.param_names_long(2) = {'phi'};
M_.param_names(3) = {'lambda'};
M_.param_names_tex(3) = {'lambda'};
M_.param_names_long(3) = {'lambda'};
M_.param_names(4) = {'beta'};
M_.param_names_tex(4) = {'beta'};
M_.param_names_long(4) = {'beta'};
M_.param_names(5) = {'sg'};
M_.param_names_tex(5) = {'sg'};
M_.param_names_long(5) = {'sg'};
M_.param_names(6) = {'rho_g'};
M_.param_names_tex(6) = {'rho\_g'};
M_.param_names_long(6) = {'rho_g'};
M_.param_names(7) = {'phi_pi'};
M_.param_names_tex(7) = {'phi\_pi'};
M_.param_names_long(7) = {'phi_pi'};
M_.param_names(8) = {'sigma_g'};
M_.param_names_tex(8) = {'sigma\_g'};
M_.param_names_long(8) = {'sigma_g'};
M_.param_names(9) = {'omega'};
M_.param_names_tex(9) = {'omega'};
M_.param_names_long(9) = {'omega'};
M_.param_names(10) = {'phi_ygap'};
M_.param_names_tex(10) = {'phi\_ygap'};
M_.param_names_long(10) = {'phi_ygap'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 5;
M_.param_nbr = 10;
M_.orig_endo_nbr = 5;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.eq_nbr = 5;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 2 7;
 0 3 8;
 0 4 0;
 0 5 0;
 1 6 0;]';
M_.nstatic = 2;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 1;
M_.ndynamic   = 3;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , 'pi' ;
  2 , 'name' , 'ygap' ;
  3 , 'name' , 'rn' ;
  4 , 'name' , 'i' ;
  5 , 'name' , 'g' ;
};
M_.mapping.pi.eqidx = [1 2 4 ];
M_.mapping.ygap.eqidx = [1 2 4 ];
M_.mapping.i.eqidx = [2 4 ];
M_.mapping.rn.eqidx = [2 3 ];
M_.mapping.g.eqidx = [3 5 ];
M_.mapping.u_g.eqidx = [5 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 2;
M_.block_structure.block(1).mfs = 2;
M_.block_structure.block(1).equation = [ 5 3];
M_.block_structure.block(1).variable = [ 5 4];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 4;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 3 4 ];
M_.block_structure.block(2).Simulation_Type = 7;
M_.block_structure.block(2).endo_nbr = 3;
M_.block_structure.block(2).mfs = 3;
M_.block_structure.block(2).equation = [ 2 4 1];
M_.block_structure.block(2).variable = [ 3 2 1];
M_.block_structure.block(2).is_linear = true;
M_.block_structure.block(2).NNZDerivatives = 10;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [1 2 3 0 0 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([1 2 1 2 3 2 3 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 1 2 2 2 3 3 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 3 6 8 ]);
M_.block_structure.variable_reordered = [ 5 4 3 2 1];
M_.block_structure.equation_reordered = [ 5 3 2 4 1];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 5 5;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 1;
 1 2;
 2 2;
 2 3;
 2 4;
 3 4;
 3 5;
 4 1;
 4 2;
 4 3;
 5 5;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 1;
 2 1;
 2 2;
];
M_.block_structure.dyn_tmp_nbr = 0;
M_.state_var = [5 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(5, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(10, 1);
M_.endo_trends = struct('deflator', cell(5, 1), 'log_deflator', cell(5, 1), 'growth_factor', cell(5, 1), 'log_growth_factor', cell(5, 1));
M_.NNZDerivatives = [16; 0; -1; ];
M_.dynamic_g1_sparse_rowval = int32([5 1 4 1 2 4 2 4 2 3 3 5 1 2 2 5 ]);
M_.dynamic_g1_sparse_colval = int32([5 6 6 7 7 7 8 8 9 9 10 10 11 11 12 16 ]);
M_.dynamic_g1_sparse_colptr = int32([1 1 1 1 1 2 4 7 9 11 13 15 16 16 16 16 17 ]);
M_.dynamic_g2_sparse_indices = int32([]);
M_.lhs = {
'pi'; 
'ygap'; 
'rn'; 
'i'; 
'g'; 
};
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 5];
M_.block_structure_stat.block(1).variable = [ 5];
M_.block_structure_stat.block(2).Simulation_Type = 1;
M_.block_structure_stat.block(2).endo_nbr = 1;
M_.block_structure_stat.block(2).mfs = 1;
M_.block_structure_stat.block(2).equation = [ 3];
M_.block_structure_stat.block(2).variable = [ 4];
M_.block_structure_stat.block(3).Simulation_Type = 6;
M_.block_structure_stat.block(3).endo_nbr = 3;
M_.block_structure_stat.block(3).mfs = 3;
M_.block_structure_stat.block(3).equation = [ 2 4 1];
M_.block_structure_stat.block(3).variable = [ 3 1 2];
M_.block_structure_stat.variable_reordered = [ 5 4 3 1 2];
M_.block_structure_stat.equation_reordered = [ 5 3 2 4 1];
M_.block_structure_stat.incidence.sparse_IM = [
 1 1;
 1 2;
 2 1;
 2 3;
 2 4;
 3 4;
 3 5;
 4 1;
 4 2;
 4 3;
 5 5;
];
M_.block_structure_stat.tmp_nbr = 0;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 2 1 2 3 2 3 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 1 2 2 2 3 3 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 3 6 8 ]);
M_.static_g1_sparse_rowval = int32([1 2 4 1 4 2 4 2 3 3 5 ]);
M_.static_g1_sparse_colval = int32([1 1 1 2 2 3 3 4 4 5 5 ]);
M_.static_g1_sparse_colptr = int32([1 4 6 8 10 12 ]);
load paramfile; 
set_param_value('gamma',gamma);
set_param_value('phi',phi);
set_param_value('lambda',lambda);
set_param_value('beta',beta);
set_param_value('sg',sg);
set_param_value('rho_g',rho_g);
set_param_value('phi_pi',phi_pi);
set_param_value('sigma_g',sigma_g);
set_param_value('omega',omega);
set_param_value('phi_ygap',phi_ygap);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(8)^2;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 20;
options_.order = 2;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save nkm_3eq.mat


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'nkm_3eq_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
