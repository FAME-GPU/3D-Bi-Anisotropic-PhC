function init_time_no_inv()
    % 初始化全局变量time_minres结构体
    global time_minres;  % 声明全局变量
    
    time_minres = struct();
    
    % 初始化minres相关字段
    time_minres.minres = 0;
    time_minres.compute_an_eig = 0;   
    time_minres.mtx_prod_vec.all = 0;
    time_minres.mtx_prod_vec.B_x = 0;
    time_minres.mtx_prod_vec.Pr_x = 0;
    time_minres.mtx_prod_vec.Qr_x = 0;
    time_minres.mtx_prod_vec.Prs_x = 0;
    time_minres.mtx_prod_vec.Qrs_x = 0;

    time_minres.sira.eig = 0;
    time_minres.sira.Choose_Ritz_Pair = 0;
      
end