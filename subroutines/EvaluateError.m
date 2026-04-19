function [res_NFSEP, res_ORGSEP] = EvaluateError(ev, ew, Eigenmode, Eigenfreq, opt, fnc, mtx_B)
    Dim = opt.comp.Dim;
    mesh_len = opt.comp.discrete.mesh_len;
    delta.x = mesh_len(1); delta.y = mesh_len(2); delta.z = mesh_len(3);
    row.x = Dim.n1; row.y = Dim.n2; row.z = Dim.n3;
    col.x = Dim.m1; col.y = Dim.m2; col.z = Dim.m3;

    C = DistOp.C(row, col, delta, opt.phys.BC);
    M = mtx_B.M_int;
    N = mtx_B.N_int;
    
    res_NFSEP  = zeros(length(Eigenfreq),1);
    res_ORGSEP = zeros(length(Eigenfreq),1);
    for j = 1:length(Eigenfreq)
        x = ev(:,j);
        lambda = ew(j);
        w = Eigenfreq(j);
        e = [Eigenmode.e1(:,j); Eigenmode.e2(:,j); Eigenmode.e3(:,j)];
        d = fnc.invN(e);

        res_NFSEP(j)  = gather( norm( fnc.Ar(x) - lambda\x )/norm(fnc.Ar(x)) );
        res_ORGSEP(j) = gather( norm( C'*(M*(C*(N*d))) - (w^2)*d ) );
    end

end