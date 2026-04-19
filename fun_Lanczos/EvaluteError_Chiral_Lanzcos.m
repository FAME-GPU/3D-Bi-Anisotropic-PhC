function [res_NFSEP, res_ORGSEP] = EvaluteError_Chiral_Lanzcos(ev, ew, Eigenmode, opt, fnc_A, fnc_B, mtx_B)
global flag_gpu

    Dim = opt.comp.Dim;
    mesh_len = opt.comp.discrete.mesh_len;
    delta.x = mesh_len(1); delta.y = mesh_len(2); delta.z = mesh_len(3);
    row.x = Dim.n1; row.y = Dim.n2; row.z = Dim.n3;
    col.x = Dim.m1; col.y = Dim.m2; col.z = Dim.m3;
    C = DistOp.C(row, col, delta, opt.phys.BC);
	
	
	M = mtx_B.M_int;
	N = mtx_B.N_int;
	Z = mtx_B.L_int;  % zeta
	X = mtx_B.J_int;  % xi
	
	fnc_mtx_B = @(x) 1i * [-Z*x(1:end/2) - M*x(end/2+1:end); N*x(1:end/2) + X * x(end/2+1:end)];	
	fnc_error_ORSEP = @(x, omega) norm( [C*x(1:end/2); C'*x(end/2+1:end)] - omega * fnc_mtx_B(x) );	
	fnc_error_NFSEP = @(x, omega) norm(fnc_A(x) - (1/omega)*fnc_B(x));
    
    res_NFSEP  = zeros(length(ew),1);
    res_ORGSEP = zeros(length(ew),1);
    if flag_gpu
        res_NFSEP = gpuArray(res_NFSEP);
        res_ORGSEP = gpuArray(res_ORGSEP);
    end
    
    for j = 1:length(ew)
        % res_NFSEP(j) = norm( fnc_A(ev(:, j)) - ew(j) * fnc_B(ev(:, j)) ) / norm(fnc.Ar(x));
        e = [Eigenmode.e1(:,j); Eigenmode.e2(:,j); Eigenmode.e3(:,j)];
        h = [Eigenmode.h1(:,j); Eigenmode.h2(:,j); Eigenmode.h3(:,j)];
        e_h = [e;h];
        e_h = e_h / norm(e_h);
        res_NFSEP(j) = fnc_error_NFSEP(ev(:, j), ew(j));
		res_ORGSEP(j) = fnc_error_ORSEP(e_h, ew(j));
        % res_ORGSEP(j) = norm( 1i * [ -1 * C' * h; C * e] - ew(j) * [ N\e + J\h; L\e + M\h] );
    end
end