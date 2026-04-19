%% compute maximum component of eigenvector

clear
clc
close all

n1 = 60;
n2 = n1;
n3 = n2;
n = n1 * n2 * n3;

%% epsilon = '2 by 2'; 
load ev_60_1.23399.mat % for 2 by 2 epsilon
vec = ev(40 * n + 1 : 41 * n);
meoa = []; %% max_ev_out_all
for i = 0.1 : 0.1 : 3.0
    str_i = num2str(i);
    load(['FCC_' str_i '.mat']);
    vec_in = vec(face_idx_in.X);
    face_idx_out.X = setdiff(1 : n, face_idx_in.X);
    vec_out = vec(face_idx_out.X);
    max_ev_out = max(abs(vec_out));
    meoa = [meoa; max_ev_out];      
end
save('max_ev_out_2by2.mat', 'meoa');
figure
semilogy(0.1 : 0.1 : 3.0, meoa, 'b*-')
% title('m_h(\rho) vs. various \rho');
xlabel('Ratio \rho')
ylabel('Maximal norm m_h(\rho)')


%% epsilon = '3 by 3'; 
load EH_SiO2_60_1.228142.mat % for 3 by 3 epsilon
vec = ev(6 * n + 1 : 7 * n);
meoa = []; %% max_ev_out_all
for i = 0.1 : 0.1 : 3.0  
    str_i = num2str(i);
    load(['FCC_' str_i '.mat']);
    vec_in = vec(edge_idx_in.Z);
    edge_idx_out.Z = setdiff(1 : n, edge_idx_in.Z);
    vec_out = vec(edge_idx_out.Z);
    max_ev_out = max(abs(vec_out));
    meoa = [meoa; max_ev_out];   
end
save('max_ev_out_3by3.mat', 'meoa');
figure
semilogy(0.1 : 0.1 : 3.0, meoa, 'b*-')
% title('m_e(\rho) vs. various \rho');
xlabel('Ratio \rho')
ylabel('Maximal norm m_e(\rho)')

