clear
clear global
clc
close all
global time_pcg
global time_FFT
global flag_gpu
global time_minres
time_FFT = 0;
addpath('subroutines')
addpath(fullfile(pwd,'InputData','BS_Gyroid'))
addpath('light_locked_Huang') 
addpath('Chiral_Media')
addpath('Jacobi_Davidson')
addpath('Public') 
addpath('fun_no_inv')
%% Load Physical and Computational settings
ew_tmp = 0; % store the last eigenvalue
init_time_no_inv(); 
flag_gpu = 1;
if flag_gpu
    gpuDevice(1);
    ew_tmp = gpuArray(ew_tmp);
end
N = 60;
gamma = 1.009; 
intep_type = 'bianisotripic'; 
model_type = 'FCC_tri'; 
filename = sprintf('%s_n%d_gamma_choose4_%g.txt', model_type, N, gamma);
fid = fopen(filename, 'w+');
fprintf("+++++++++++++++ gamma = %.8d, begin +++++++++++++\n", gamma);
for radius_rate = 1 %.217:0.001:1.3
[ opt.phys ] = Physical_settings_tri(model_type, gamma, radius_rate); % HEX
[ opt.comp ] = Computational_settings(opt.phys, N);
opt.comp.eigensolver.eigen_wanted = 5;
if 0
    Plot_Material_Isosurf(opt);
end
%% Construct discrete material parameter matrices
[ mtx_B ] = Material_Matrix_Generate_Mult_Point( opt, intep_type,radius_rate ); % material related matrices generate

    if 0 % load grid (radius)
        fprintf("+++++++++++++++ radius_rate = %.8d +++++++++++++\n", radius_rate);
        Standard_inner_idx = mtx_B.Standard_inner_idx;
        Bodycenter_inner_idx = mtx_B.Bodycenter_inner_idx;
        filename = sprintf('../手性/FCC_Sphere_radius=0.2/radius_data/radius_rate=%.2f.mat', radius_rate);
        save(filename, 'Standard_inner_idx', 'Bodycenter_inner_idx');
    end        
end
result = cell(1,size(opt.phys.lattice.wave_vec_array,2));
Freq_mtx = zeros(opt.comp.eigensolver.eigen_wanted, size(opt.phys.lattice.wave_vec_array,2));
if flag_gpu
    Freq_mtx = gpuArray(Freq_mtx);
end
i = 51;
result{i}.opt = opt;
result{i}.opt.phys.lattice.wave_vec_idx = i;
wave_vec = opt.phys.lattice.wave_vec_array(:,i);
if norm(wave_vec) == 0
    wave_vec = wave_vec + (1e-4);
end
lattice_vec = opt.phys.lattice.lattice_vector;
opt.phys.BC.x.theta = dot( wave_vec, lattice_vec(:,1) );
opt.phys.BC.y.theta = dot( wave_vec, lattice_vec(:,2) );
opt.phys.BC.z.theta = dot( wave_vec, lattice_vec(:,3) );
%% Construct SVD matrices, interploated material matrices, and function handles
[ mtx_A ] = SVD_Matrix_Generate( opt ); % SVD matrices generate
[ mtx_B ] = Interpolation_Material_Matrix_Mult_Point_vertex(opt, mtx_B, intep_type); % Interpolate material matrices with positive definiteness

if opt.comp.flag_gpu
    [ mtx_A, mtx_B ] = Convert2gpuarray(mtx_A,mtx_B);
end
[ fnc   ] = Funhand_Generate( opt.comp.Dim, mtx_A, mtx_B, opt.comp.linearsolver, opt.comp.NFSEP_type );

%% Solving NFSEP
fprintf('Start solving NFSEP for %d eigenvalues\n', opt.comp.eigensolver.eigen_wanted)

%% double-curl && uncoupled && epsilon is singular
length_sigma_r = length(mtx_A.Lambda_curl.Sigma.r);
mtx_A.Lambda_curl.Sigma.inv_r = 1./mtx_A.Lambda_curl.Sigma.r;
% 一些准备函数
fncA = @(x) 1i*[mtx_A.Lambda_curl.Sigma.inv_r .* x(end/2+1:end); ...
    -mtx_A.Lambda_curl.Sigma.inv_r .* x(1:end/2)];
fncB = @(y)fncB_Generate_no_inv(fnc,mtx_B,y);

kind = 'SIRA';  % Lanczos  SIRA
if  strcmp(kind, 'Lanczos')
    %% Lanczos
    op_A = @(x)fncA(x, mtx_A);
    op_B = @(x)fncB(x, fnc, mtx_B);
    op_invB = @(x)fncinvB(x, op_B, 1e-12, 1e4);

    v0 = ones(2 * length(mtx_A.Lambda_curl.Sigma.r),1);
    vB = fncB(v0, fnc, mtx_B);
    v0 = v0/sqrt(abs(dot(v0,vB)));
    [ev, ew, current_iter] = Lanczos(op_A, op_B, op_invB, 4 * N^3, opt.comp.eigensolver.eigen_wanted, opt.comp.eigensolver.tol, opt.comp.eigensolver.itmax, opt.comp.eigensolver.p, v0, opt.comp.flag_gpu);
    iter_records(i) = current_iter;
    ew = 1./ew;
    Freq_mtx(:,i) = ew;
    [ result{i}.Eigenmode.e1, result{i}.Eigenmode.e2, result{i}.Eigenmode.e3, ...
    result{i}.Eigenmode.h1, result{i}.Eigenmode.h2, result{i}.Eigenmode.h3] ...
     = Eigen_Restoration_Chiral_no_inv(ev, opt.comp.Dim, mtx_B, fnc.Qr, fnc.Pr);
    [result{i}.res_NFSEP, result{i}.res_ORGSEP] = EvaluteError_Chiral_no_inv(ev, ew, result{i}.Eigenmode, opt, op_A, op_B, mtx_B);
    fprintf('Done: %d/%d.\n', i, size(opt.phys.lattice.wave_vec_array,2));
elseif strcmp(kind, 'SIRA')
    
    initial_V = rand(2 * length_sigma_r, 1, 'double', 'gpuArray') + 1i*rand(2 * length_sigma_r, 1, 'double', 'gpuArray');
    initial_V = initial_V / norm(initial_V);
    if flag_gpu
        initial_V = gpuArray(initial_V);
    end
    
    stop_tolerance = 1.0e-8; 
    target_type    = 'CL';
    
    sigma = 0;
    ncv = opt.comp.eigensolver.eigen_wanted;
    no_restart     = 50; %35;
    mtxdim         = 2*length(mtx_A.Lambda_curl.Sigma.r);  % n=6;  mtxdim         = 32 * n
    CorEq          = 'SIRA'; %'SIRA'; %'JD'; %'SIRA';
    LSinfo.solver  = 'minres';%'gmres'; %'bicgstabl'; %'minres';%'bicgstabl'; %'bicgstabl'; %'minres';
    precond = 'no';
    LSinfo.precond = precond; %'no'; %'yes';
    LSinfo.tol = 1e-2;
    LSinfo.maxit   = 30000;
    [conv_ew, conv_ev] = GEP_AB_Herm_JDSIRA_Driver4 (fncB, fncA, mtxdim, no_restart, ...
        ncv, stop_tolerance, initial_V, sigma, fid, target_type, CorEq, LSinfo, [], [], opt, mtx_B, fnc, wave_vec, N, gamma, model_type);
    ew = conv_ew;
end
if ( abs(imag(ew(1))) > 1e-3 || abs(ew_tmp-ew(1))/abs(ew(1)) >1e-2)
    fprintf("======================== gamma = %.8d, ew: %d+%di ============\n", gamma, real(ew(1)), imag(ew(1)));
end
ew_tmp = ew(1);
time_pcg.invtildeN.Phis(1) = [];
time_pcg.invtildeN.Phi(1) = [];
time_pcg.invtildeM.Psis(1) = [];
time_pcg.invtildeM.Psi(1) = [];
result{i}.EventTime.time_pcg = time_pcg;

if strcmp(kind, 'Lanczos')
    Freq_mtx(:,i) = abs(real(ew));
elseif strcmp(kind, 'SIRA')
    Freq_mtx(:,i) = ew;
end

fprintf('Done: %d/%d.\n', i, size(opt.phys.lattice.wave_vec_array,2));
fclose(fid);
