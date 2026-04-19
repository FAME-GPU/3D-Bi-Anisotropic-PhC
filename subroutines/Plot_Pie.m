function Plot_Pie(result)
    % 
    time_Psis = EventTime.Psi('t');  time_Psis.reshape = 0; time_Psis = rmfield(time_Psis,{'rs1','rs2','rs3'});
    time_Psi  = EventTime.Psi('nt'); time_Psi.reshape  = 0; time_Psi  = rmfield(time_Psi ,{'rs1','rs2','rs3'});
    time_Phis = EventTime.Phi('t');  time_Phis.reshape = 0; time_Phis = rmfield(time_Phis,{'rs1','rs2','rs3'});
    time_Phi  = EventTime.Phi('nt'); time_Phi.reshape  = 0; time_Phi  = rmfield(time_Phi ,{'rs1','rs2','rs3'});
    
    count = [0 0];
    % 
    for i = 1:length(result)
        tmp_M = result{i}.EventTime.time_pcg.invtildeM;
        tmp_N = result{i}.EventTime.time_pcg.invtildeN;
        %%
        count(1) = count(1) + sum(tmp_M.iter);
        count(2) = count(2) + sum(tmp_N.iter);

        time_Psis.D    = time_Psis.D    + sum(vertcat(tmp_M.Psis.D   ));
        time_Psis.reshape = time_Psis.reshape  + sum(vertcat(tmp_M.Psis.rs1 ))...
                                               + sum(vertcat(tmp_M.Psis.rs2 ))...
                                               + sum(vertcat(tmp_M.Psis.rs3 ));
        time_Psis.P1   = time_Psis.P1   + sum(vertcat(tmp_M.Psis.P1  ));
        time_Psis.P2   = time_Psis.P2   + sum(vertcat(tmp_M.Psis.P2  ));
        time_Psis.P3   = time_Psis.P3   + sum(vertcat(tmp_M.Psis.P3  ));
        time_Psis.fft1 = time_Psis.fft1 + sum(vertcat(tmp_M.Psis.fft1));
        time_Psis.fft2 = time_Psis.fft2 + sum(vertcat(tmp_M.Psis.fft2));
        time_Psis.fft3 = time_Psis.fft3 + sum(vertcat(tmp_M.Psis.fft3));
        time_Psis.R1   = time_Psis.R1   + sum(vertcat(tmp_M.Psis.R1  ));
        time_Psis.R2   = time_Psis.R2   + sum(vertcat(tmp_M.Psis.R2  ));
        time_Psis.R3   = time_Psis.R3   + sum(vertcat(tmp_M.Psis.R3  ));
        time_Psis.cat  = time_Psis.cat  + sum(vertcat(tmp_M.Psis.cat ));
        time_Psis.S    = time_Psis.S    + sum(vertcat(tmp_M.Psis.S   ));
        %%
        time_Psi.D     = time_Psi.D     + sum(vertcat(tmp_M.Psi.D   ));
        time_Psi.reshape = time_Psi.reshape  + sum(vertcat(tmp_M.Psi.rs1 ))...
                                             + sum(vertcat(tmp_M.Psi.rs2 ))...
                                             + sum(vertcat(tmp_M.Psi.rs3 ));
        time_Psi.P1    = time_Psi.P1    + sum(vertcat(tmp_M.Psi.P1  ));
        time_Psi.P2    = time_Psi.P2    + sum(vertcat(tmp_M.Psi.P2  ));
        time_Psi.P3    = time_Psi.P3    + sum(vertcat(tmp_M.Psi.P3  ));
        time_Psi.ifft1 = time_Psi.ifft1 + sum(vertcat(tmp_M.Psi.ifft1));
        time_Psi.ifft2 = time_Psi.ifft2 + sum(vertcat(tmp_M.Psi.ifft2));
        time_Psi.ifft3 = time_Psi.ifft3 + sum(vertcat(tmp_M.Psi.ifft3));
        time_Psi.R1    = time_Psi.R1    + sum(vertcat(tmp_M.Psi.R1  ));
        time_Psi.R2    = time_Psi.R2    + sum(vertcat(tmp_M.Psi.R2  ));
        time_Psi.R3    = time_Psi.R3    + sum(vertcat(tmp_M.Psi.R3  ));
        time_Psi.cat   = time_Psi.cat   + sum(vertcat(tmp_M.Psi.cat ));
        time_Psi.S     = time_Psi.S     + sum(vertcat(tmp_M.Psi.S   ));
        %%
        time_Phis.D    = time_Phis.D    + sum(vertcat(tmp_N.Phis.D   ));
        time_Phis.reshape = time_Phis.reshape  + sum(vertcat(tmp_N.Phis.rs1 ))...
                                               + sum(vertcat(tmp_N.Phis.rs2 ))...
                                               + sum(vertcat(tmp_N.Phis.rs3 ));
        time_Phis.P1   = time_Phis.P1   + sum(vertcat(tmp_N.Phis.P1  ));
        time_Phis.P2   = time_Phis.P2   + sum(vertcat(tmp_N.Phis.P2  ));
        time_Phis.P3   = time_Phis.P3   + sum(vertcat(tmp_N.Phis.P3  ));
        time_Phis.fft1 = time_Phis.fft1 + sum(vertcat(tmp_N.Phis.fft1));
        time_Phis.fft2 = time_Phis.fft2 + sum(vertcat(tmp_N.Phis.fft2));
        time_Phis.fft3 = time_Phis.fft3 + sum(vertcat(tmp_N.Phis.fft3));
        time_Phis.R1   = time_Phis.R1   + sum(vertcat(tmp_N.Phis.R1  ));
        time_Phis.R2   = time_Phis.R2   + sum(vertcat(tmp_N.Phis.R2  ));
        time_Phis.R3   = time_Phis.R3   + sum(vertcat(tmp_N.Phis.R3  ));
        time_Phis.cat  = time_Phis.cat  + sum(vertcat(tmp_N.Phis.cat ));
        time_Phis.S    = time_Phis.S    + sum(vertcat(tmp_N.Phis.S   ));
        %%
        time_Phi.D     = time_Phi.D     + sum(vertcat(tmp_N.Phi.D   ));
        time_Phi.reshape = time_Phi.reshape  + sum(vertcat(tmp_N.Phi.rs1 ))...
                                             + sum(vertcat(tmp_N.Phi.rs2 ))...
                                             + sum(vertcat(tmp_N.Phi.rs3 ));
        time_Phi.P1    = time_Phi.P1    + sum(vertcat(tmp_N.Phi.P1  ));
        time_Phi.P2    = time_Phi.P2    + sum(vertcat(tmp_N.Phi.P2  ));
        time_Phi.P3    = time_Phi.P3    + sum(vertcat(tmp_N.Phi.P3  ));
        time_Phi.ifft1 = time_Phi.ifft1 + sum(vertcat(tmp_N.Phi.ifft1));
        time_Phi.ifft2 = time_Phi.ifft2 + sum(vertcat(tmp_N.Phi.ifft2));
        time_Phi.ifft3 = time_Phi.ifft3 + sum(vertcat(tmp_N.Phi.ifft3));
        time_Phi.R1    = time_Phi.R1    + sum(vertcat(tmp_N.Phi.R1  ));
        time_Phi.R2    = time_Phi.R2    + sum(vertcat(tmp_N.Phi.R2  ));
        time_Phi.R3    = time_Phi.R3    + sum(vertcat(tmp_N.Phi.R3  ));
        time_Phi.cat   = time_Phi.cat   + sum(vertcat(tmp_N.Phi.cat ));
        time_Phi.S     = time_Phi.S     + sum(vertcat(tmp_N.Phi.S   ));
       
    end
    figure(1)
    tiledlayout(1,2,'TileSpacing','compact');
    ax1 = nexttile;
    X = cell2mat(struct2cell(time_Psi)); t1 = sum(X); X = X / t1; 
    pie(ax1,X)
    title(['\psi (avg. time =',num2str(t1/count(1),'%.2e'),' sec.)'])
    ax2 = nexttile;
    X = cell2mat(struct2cell(time_Phi)); t2 = sum(X); X = X / t2;
    pie(ax2,X,fieldnames(time_Phi))
    title(['\phi (avg. time =',num2str(t2/count(2),'%.2e'),' sec.)'])
    lgd = legend(fieldnames(time_Phi));
    lgd.Layout.Tile = 'east';
    
    figure(2)
    tiledlayout(1,2,'TileSpacing','compact');
    ax3 = nexttile;
    X = cell2mat(struct2cell(time_Psis)); t3 = sum(X); X = X / t3;
    pie(ax3,X,fieldnames(time_Psis))
    title(['\psi^* (avg. time =',num2str(t3/count(1),'%.2e'),' sec.)'])
    ax4 = nexttile;
    X = cell2mat(struct2cell(time_Phis)); t4 = sum(X); X = X / t4;
    pie(ax4,X,fieldnames(time_Phis))
    title(['\phi^* (avg. time =',num2str(t4/count(2),'%.2e'),' sec.)'])

    lgd = legend(fieldnames(time_Psis));
    lgd.Layout.Tile = 'east';
end