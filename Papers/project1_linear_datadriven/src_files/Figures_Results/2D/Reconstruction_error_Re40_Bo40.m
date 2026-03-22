clc; clear; close all

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
POD_path  = fullfile(scriptDir, '..', '..', 'Codes', 'POD', 'Results');
POD_path  = char(java.io.File(POD_path).getCanonicalPath());

make_pod = @(obs, s, e, DR, VR, Re, Bo, var) ...
    fullfile(POD_path, sprintf( ...
        'Results_POD_%s_2Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
        obs, s, e, DR, VR, Re, Bo, var));

DMD_path  = fullfile(scriptDir, '..', '..', 'Codes', 'DMD', 'Results');
DMD_path  = char(java.io.File(DMD_path).getCanonicalPath());

make_dmd = @(obs, s, e, DR, VR, Re, Bo, var) ...
    fullfile(DMD_path, sprintf( ...
        'Results_DMD_%s_2Dnew_%d_%d_DR%d_VR%d_Re_%d_Bo_%d_%s.h5', ...
        obs, s, e, DR, VR, Re, Bo, var));

%% Case setup
start_num = 1;
end_num   = 400;
DR        = 10;
VR        = 10;
Re        = 40;
Bo        = 40;
nmax      = 50;

grp_sigma = '/Raw/Sigma';
grp_re    = '/Reconstruction_error';
grp_rv    = '/Standard_Decomposition/Recons_vec';
grp_svd   = '/Standard_Decomposition/SVD_S';

%% POD mean1 and mean0 sigma
EG_m1   = h5read(make_pod('EG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_sigma);
EV_m1   = h5read(make_pod('EUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_sigma);
LG_m1   = h5read(make_pod('LG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_sigma);
LV_m1   = h5read(make_pod('LUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_sigma);
LEG_m1  = h5read(make_pod('LEGi',  start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_sigma);
LEV_m1  = h5read(make_pod('LEUVi', start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_sigma);
LP_d1   = h5read(make_pod('LP',    start_num,end_num,DR,VR,Re,Bo,'der1'),  grp_sigma);

EG_m0   = h5read(make_pod('EG',    start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_sigma);
EV_m0   = h5read(make_pod('EUV',   start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_sigma);
LG_m0   = h5read(make_pod('LG',    start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_sigma);
LV_m0   = h5read(make_pod('LUV',   start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_sigma);
LEG_m0  = h5read(make_pod('LEGi',  start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_sigma);
LEV_m0  = h5read(make_pod('LEUVi', start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_sigma);

%% POD reconstruction errors
re_EG   = h5read(make_pod('EG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_re);
re_EV   = h5read(make_pod('EUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_re);
re_LG   = h5read(make_pod('LG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_re);
re_LV   = h5read(make_pod('LUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_re);
re_LEG  = h5read(make_pod('LEGi',  start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_re);
re_LEV  = h5read(make_pod('LEUVi', start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_re);
re_LP   = h5read(make_pod('LP',    start_num,end_num,DR,VR,Re,Bo,'der1'),  grp_re);

%% DMD recons vec and SVD_S
rv_EG   = h5read(make_dmd('EG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_rv);
rv_EV   = h5read(make_dmd('EUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_rv);
rv_LG   = h5read(make_dmd('LG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_rv);
rv_LV   = h5read(make_dmd('LUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_rv);
rv_LEG  = h5read(make_dmd('LEGi',  start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_rv);
rv_LEV  = h5read(make_dmd('LEUVi', start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_rv);
rv_LP   = h5read(make_dmd('LP',    start_num,end_num,DR,VR,Re,Bo,'der1'),  grp_rv);

ss0_EG  = h5read(make_dmd('EG',    start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_svd);
ss0_EV  = h5read(make_dmd('EUV',   start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_svd);
ss0_LG  = h5read(make_dmd('LG',    start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_svd);
ss0_LV  = h5read(make_dmd('LUV',   start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_svd);
ss0_LEG = h5read(make_dmd('LEGi',  start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_svd);
ss0_LEV = h5read(make_dmd('LEUVi', start_num,end_num,DR,VR,Re,Bo,'mean0'), grp_svd);

ss1_EG  = h5read(make_dmd('EG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_svd);
ss1_EV  = h5read(make_dmd('EUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_svd);
ss1_LG  = h5read(make_dmd('LG',    start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_svd);
ss1_LV  = h5read(make_dmd('LUV',   start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_svd);
ss1_LEG = h5read(make_dmd('LEGi',  start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_svd);
ss1_LEV = h5read(make_dmd('LEUVi', start_num,end_num,DR,VR,Re,Bo,'mean1'), grp_svd);

%% Fair reconstruction errors
fe_EG  = re_EG  * (sum(abs(EG_m1))  / sum(abs(EG_m0)));
fe_EV  = re_EV  * (sum(abs(EV_m1))  / sum(abs(EV_m0)));
fe_LG  = re_LG  * (sum(abs(LG_m1))  / sum(abs(LG_m0)));
fe_LV  = re_LV  * (sum(abs(LV_m1))  / sum(abs(LV_m0)));
fe_LEG = re_LEG * (sum(abs(LEG_m1)) / sum(abs(LEG_m0)));
fe_LEV = re_LEV * (sum(abs(LEV_m1)) / sum(abs(LEV_m0)));
fe_LP  = re_LP;

dfe_EG  = dmd_err_vec(rv_EG,  nmax) * (sum(diag(ss1_EG).^2)  / sum(diag(ss0_EG).^2));
dfe_EV  = dmd_err_vec(rv_EV,  nmax) * (sum(diag(ss1_EV).^2)  / sum(diag(ss0_EV).^2));
dfe_LG  = dmd_err_vec(rv_LG,  nmax) * (sum(diag(ss1_LG).^2)  / sum(diag(ss0_LG).^2));
dfe_LV  = dmd_err_vec(rv_LV,  nmax) * (sum(diag(ss1_LV).^2)  / sum(diag(ss0_LV).^2));
dfe_LEG = dmd_err_vec(rv_LEG, nmax) * (sum(diag(ss1_LEG).^2) / sum(diag(ss0_LEG).^2));
dfe_LEV = dmd_err_vec(rv_LEV, nmax) * (sum(diag(ss1_LEV).^2) / sum(diag(ss0_LEV).^2));
dfe_LP  = dmd_err_vec(rv_LP,  nmax);

%% Console summary at r = 10
fprintf('\nRe40 Bo40 fair GVF reconstruction errors at r = 10\n');
fprintf('Eulerian  POD = %.4f%%, DMD = %.4f%%\n', fe_EG(10),  dfe_EG(10));
fprintf('Lagrangian POD = %.4f%%, DMD = %.4f%%\n', fe_LG(10),  dfe_LG(10));
fprintf('Moving-frame POD = %.4f%%, DMD = %.4f%%\n', fe_LEG(10), dfe_LEG(10));

fprintf('\nRe40 Bo40 velocity reconstruction errors at r = 10\n');
fprintf('Eulerian velocity POD = %.4f%%, DMD = %.4f%%\n', fe_EV(10), dfe_EV(10));
fprintf('Lagrangian velocity POD = %.4f%%, DMD = %.4f%%\n', fe_LV(10), dfe_LV(10));
fprintf('Moving-frame velocity POD = %.4f%%, DMD = %.4f%%\n', fe_LEV(10), dfe_LEV(10));
fprintf('Lagrangian position POD = %.4f%%, DMD = %.4f%%\n', fe_LP(10), dfe_LP(10));

%% Export arrays for later manuscript analysis
export_data = struct();
export_data.rank = (1:nmax)';

export_data.pod.gvf.re40.eulerian = fe_EG(1:nmax);
export_data.pod.gvf.re40.lagrangian = fe_LG(1:nmax);
export_data.pod.gvf.re40.moving_frame = fe_LEG(1:nmax);
export_data.pod.velocity.re40.eulerian = fe_EV(1:nmax);
export_data.pod.velocity.re40.lagrangian = fe_LV(1:nmax);
export_data.pod.velocity.re40.moving_frame = fe_LEV(1:nmax);
export_data.pod.velocity.re40.lagrangian_position = fe_LP(1:nmax);
export_data.dmd.gvf.re40.eulerian = dfe_EG(1:nmax);
export_data.dmd.gvf.re40.lagrangian = dfe_LG(1:nmax);
export_data.dmd.gvf.re40.moving_frame = dfe_LEG(1:nmax);
export_data.dmd.velocity.re40.eulerian = dfe_EV(1:nmax);
export_data.dmd.velocity.re40.lagrangian = dfe_LV(1:nmax);
export_data.dmd.velocity.re40.moving_frame = dfe_LEV(1:nmax);
export_data.dmd.velocity.re40.lagrangian_position = dfe_LP(1:nmax);

metadata = struct();
metadata.nmax = nmax;
metadata.note = 'Arrays exported from Reconstruction_error_Re40_Bo40.m. These are fair reconstruction-error profiles for Re40 Bo40 only.';
metadata.re40 = struct('re', Re, 'bo', Bo, 'start_num', start_num, 'end_num', end_num);

save_reconstruction_error_export(export_data, metadata, 'reconstruction_error_analysis_Re40_Bo40_2D.mat');

%% Local function
function err_vec = dmd_err_vec(rv, nmax)
    ranks  = double(rv(2,:));
    errors = double(rv(4,:));
    [r_sorted, si] = sort(ranks);
    e_sorted = errors(si);
    err_vec = interp1(r_sorted, e_sorted, (1:nmax)', 'linear', 'extrap');
end