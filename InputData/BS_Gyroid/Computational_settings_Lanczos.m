function comp = Computational_settings_Lanczos(physical, N)
global flag_PlotMaterial gpu flag_gpu flag_FFT

flag_PlotMaterial = 1;

flag_FFT = 1; % 0: DFT/DCT/DST; 1: pure FFT; % original value: 0

comp.flag_gpu = flag_gpu;
if flag_gpu
    %comp.gpudev = gpuDevice(1);
    %gpu = comp.gpudev;
    %reset(comp.gpudev);
end
comp.flag_numerical_validation = 0;
%% NFSEP type
% comp.NFSEP_type: Solve Ar*x = λ*x, where
%   0: Ar:=M\Σr\N\Σr, x:=Σr*er, λ:=w^{-2}
%   1: Ar:=[ 0, Σr]^{-1/2}[-N, 0]^{-1}[ 0, Σr]^{-1/2}  x:=[ 0, Σr]^{1/2}[hr]  λ:=-iw^{-1}
%          [Σr,  0]       [ 0, M]     [Σr,  0]       ,    [Σr,  0][er]      , 

comp.NFSEP_type = 0;
%% Eigensolver settings
% outer iteration settings
comp.eigensolver.itmax        = 10000; % 10000
comp.eigensolver.tol          = 1e-10;
comp.eigensolver.isreal       = 0;
comp.eigensolver.issym        = 0;
comp.eigensolver.eigen_wanted = 10;
comp.eigensolver.target       = 'lm';
comp.eigensolver.p            = 2 * comp.eigensolver.eigen_wanted; % 30
% inner iteration settings
%comp.linearsolver.solver = 'pcg';
comp.linearsolver.itmax  = 1000;
comp.linearsolver.tol    = 1e-12;
comp.linearsolver.x0_fnc = @(n) zeros(n,1);
%% Discretization settings
Yee_cell_num = [N,N,N]; % origin: 4,4,4
comp.discrete.mesh_len = physical.lattice.constant ./ Yee_cell_num;
comp.discrete.grid_num = (physical.geometry.supercell_num).*Yee_cell_num;
%% Parameter setting
n1 = comp.discrete.grid_num(1);
n2 = comp.discrete.grid_num(2);
n3 = comp.discrete.grid_num(3);
m1 = n1 - strcmp(physical.BC.x.str,'PEC');
m2 = n2 - strcmp(physical.BC.y.str,'PEC');
m3 = n3 - strcmp(physical.BC.z.str,'PEC');

comp.Dim.n1 = n1;
comp.Dim.n2 = n2;
comp.Dim.n3 = n3;
comp.Dim.m1 = m1;
comp.Dim.m2 = m2;
comp.Dim.m3 = m3;

comp.Dim.rmc_1D.x = n1 - m1;
comp.Dim.rmc_1D.y = n2 - m2;
comp.Dim.rmc_1D.z = n3 - m3;

comp.Dim.E1 = [n1,m2,m3]; 
comp.Dim.E2 = [m1,n2,m3]; 
comp.Dim.E3 = [m1,m2,n3];
comp.Dim.H1 = [m1,n2,n3]; 
comp.Dim.H2 = [n1,m2,n3]; 
comp.Dim.H3 = [n1,n2,m3];

comp.Dim.Ne1 = n1 * m2 * m3;
comp.Dim.Ne2 = m1 * n2 * m3;
comp.Dim.Ne3 = m1 * m2 * n3;

comp.Dim.Nh1 = m1 * n2 * n3;
comp.Dim.Nh2 = n1 * m2 * n3;
comp.Dim.Nh3 = n1 * n2 * m3;

comp.Dim.Ne = comp.Dim.Ne1 + comp.Dim.Ne2 + comp.Dim.Ne3;
comp.Dim.Nh = comp.Dim.Nh1 + comp.Dim.Nh2 + comp.Dim.Nh3;

comp.Dim.N  = n1 * n2 * n3;
comp.Dim.Nq = m1 * m2 * m3;
    
comp.Dim.Nhat1 = m1 * comp.Dim.rmc_1D.y * comp.Dim.rmc_1D.z;
comp.Dim.Nhat2 = comp.Dim.rmc_1D.x * m2 * comp.Dim.rmc_1D.z;
comp.Dim.Nhat3 = comp.Dim.rmc_1D.x * comp.Dim.rmc_1D.y * m3;

comp.Dim.Ntil1 = comp.Dim.rmc_1D.x * m2 * m3;
comp.Dim.Ntil2 = m1 * comp.Dim.rmc_1D.y * m3;
comp.Dim.Ntil3 = m1 * m2 * comp.Dim.rmc_1D.z;

comp.Dim.Ntil = comp.Dim.Ntil1 + comp.Dim.Ntil2 + comp.Dim.Ntil3;
comp.Dim.Nhat = comp.Dim.Nhat1 + comp.Dim.Nhat2 + comp.Dim.Nhat3;
    
comp.Dim.Nr = 2*comp.Dim.Nq + comp.Dim.Ntil;
comp.Dim.Ar = (comp.NFSEP_type + 1) * comp.Dim.Nr;

comp.Dim.Phi.SIZE = [comp.Dim.E1; comp.Dim.E2; comp.Dim.E3];
comp.Dim.Psi.SIZE = [comp.Dim.H1; comp.Dim.H2; comp.Dim.H3];

comp.Dim.FFT_SIZE.x = n1*(1+strcmp(physical.BC.x.str,'PEC'));
comp.Dim.FFT_SIZE.y = n2*(1+strcmp(physical.BC.y.str,'PEC'));
comp.Dim.FFT_SIZE.z = n3*(1+strcmp(physical.BC.z.str,'PEC'));
comp.Dim.Phi.FFT_SIZE = comp.Dim.FFT_SIZE;
comp.Dim.Psi.FFT_SIZE = comp.Dim.FFT_SIZE;
