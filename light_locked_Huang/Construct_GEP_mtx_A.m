
 function [GEP_mtx_A, mtx_C] = Construct_GEP_mtx_A(n_1,n_2,n_3,vec_k,period_vec1,period_vec2,period_vec3,delta)
% 

     n_12  =  n_1 * n_2;
     n_123 = n_12 * n_3;
     
     [i_row_1, j_col_1, k_val_1] = Gamma_1(n_1,n_2,n_3,period_vec1,vec_k);
     [i_row_2, j_col_2, k_val_2] = Gamma_2(n_1,n_2,n_3,period_vec1,period_vec2,vec_k);
     [i_row_3, j_col_3, k_val_3] = Gamma_3(n_1,n_2,n_3,period_vec1,period_vec2,period_vec3,vec_k);
     
     len_1 = length(i_row_1);
     len_2 = length(i_row_2);
     len_3 = length(i_row_3);
     
     len   = 8 * ( len_1  + len_2 + len_3 );
     
     i_row  = zeros(len,1);
     j_col  = zeros(len,1);
     k_val  = zeros(len,1);
     
     k_val_1  = k_val_1 / delta(1);
     k_val_2  = k_val_2 / delta(2);
     k_val_3  = k_val_3 / delta(3);
     
%
% Setup C(C_1, C_2, C_3)
%
% ====== First block rows

% - G_3 matrix
     idx_sta                    =     1; 
     idx_end                    = len_3;  
     
     i_row( idx_sta : idx_end , 1 ) =         i_row_3;
     j_col( idx_sta : idx_end , 1 ) = n_123 + j_col_3;
     k_val( idx_sta : idx_end , 1 ) =       - k_val_3;

%   G_2 matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row(idx_sta : idx_end , 1 ) =             i_row_2;
     j_col(idx_sta : idx_end , 1 ) = 2 * n_123 + j_col_2;
     k_val(idx_sta : idx_end , 1 ) =             k_val_2;

% ======= Second block rows 

%   G_3 matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     
     i_row( idx_sta : idx_end , 1 ) = n_123 + i_row_3;
     j_col( idx_sta : idx_end , 1 ) =         j_col_3;
     k_val( idx_sta : idx_end , 1 ) =         k_val_3;
     
% - G_1 matrix         
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) =     n_123 + i_row_1;
     j_col( idx_sta : idx_end , 1 ) = 2 * n_123 + j_col_1;
     k_val( idx_sta : idx_end , 1 ) =           - k_val_1;

% ===== Third block rows

% - G_2 matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row( idx_sta : idx_end , 1 ) = 2 * n_123 + i_row_2;
     j_col( idx_sta : idx_end , 1 ) =             j_col_2;
     k_val( idx_sta : idx_end , 1 ) =           - k_val_2;
     
%   G_1 matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = 2 * n_123 + i_row_1;
     j_col( idx_sta : idx_end , 1 ) =     n_123 + j_col_1;
     k_val( idx_sta : idx_end , 1 ) =             k_val_1;
     
%      mtx_C            = sparse(i_row, j_col, k_val);
% 
%
%
%
% Setup C(-C_1', -C_2', C_3)
%

%
% ====== First block rows

% - G_3 matrix 
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     dimC                   = 3 * n_123;
     
     i_row( idx_sta : idx_end , 1 ) = dimC +         i_row_3;
     j_col( idx_sta : idx_end , 1 ) = dimC + n_123 + j_col_3;
     k_val( idx_sta : idx_end , 1 ) =       - k_val_3;

% - G_2' matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row(idx_sta : idx_end , 1 ) = dimC +             j_col_2;
     j_col(idx_sta : idx_end , 1 ) = dimC + 2 * n_123 + i_row_2;
     k_val(idx_sta : idx_end , 1 ) =             -conj(k_val_2);

% ======= Second block rows 

%   G_3 matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     
     i_row( idx_sta : idx_end , 1 ) = dimC + n_123 + i_row_3;
     j_col( idx_sta : idx_end , 1 ) = dimC +         j_col_3;
     k_val( idx_sta : idx_end , 1 ) =         k_val_3;
     
% G_1' matrix         
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = dimC +     n_123 + j_col_1;
     j_col( idx_sta : idx_end , 1 ) = dimC + 2 * n_123 + i_row_1;
     k_val( idx_sta : idx_end , 1 ) =           conj(k_val_1);

% ===== Third block rows

% G_2' matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row( idx_sta : idx_end , 1 ) = dimC + 2 * n_123 + j_col_2;
     j_col( idx_sta : idx_end , 1 ) = dimC +             i_row_2;
     k_val( idx_sta : idx_end , 1 ) =           conj(k_val_2);
     
% - G_1' matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = dimC + 2 * n_123 + j_col_1;
     j_col( idx_sta : idx_end , 1 ) = dimC +     n_123 + i_row_1;
     k_val( idx_sta : idx_end , 1 ) =            - conj(k_val_1);
     
% 
%
%
%
% Setup C(-C_1', C_2, -C_3')
%
%
% ====== First block rows

% G_3' matrix 
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     dimC2                  = 2 * dimC;
     
     i_row( idx_sta : idx_end , 1 ) = dimC2 +         j_col_3;
     j_col( idx_sta : idx_end , 1 ) = dimC2 + n_123 + i_row_3;
     k_val( idx_sta : idx_end , 1 ) =       conj(k_val_3);

%   G_2 matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row(idx_sta : idx_end , 1 ) = dimC2 +             i_row_2;
     j_col(idx_sta : idx_end , 1 ) = dimC2 + 2 * n_123 + j_col_2;
     k_val(idx_sta : idx_end , 1 ) =             k_val_2;

% ======= Second block rows 

% - G_3' matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     
     i_row( idx_sta : idx_end , 1 ) = dimC2 + n_123 + j_col_3;
     j_col( idx_sta : idx_end , 1 ) = dimC2 +         i_row_3;
     k_val( idx_sta : idx_end , 1 ) =         - conj(k_val_3);
     
% G_1' matrix         
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = dimC2 +     n_123 + j_col_1;
     j_col( idx_sta : idx_end , 1 ) = dimC2 + 2 * n_123 + i_row_1;
     k_val( idx_sta : idx_end , 1 ) =           conj(k_val_1);

% ===== Third block rows

% - G_2 matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row( idx_sta : idx_end , 1 ) = dimC2 + 2 * n_123 + i_row_2;
     j_col( idx_sta : idx_end , 1 ) = dimC2 +             j_col_2;
     k_val( idx_sta : idx_end , 1 ) =           - k_val_2;
     
% - G_1' matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = dimC2 + 2 * n_123 + j_col_1;
     j_col( idx_sta : idx_end , 1 ) = dimC2 +     n_123 + i_row_1;
     k_val( idx_sta : idx_end , 1 ) =             - conj(k_val_1);

% 
%
%
%
% Setup C(C_1, -C_2', -C_3')
%
%
% ====== First block rows

% G_3' matrix 
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     dimC3                  = 3 * dimC;
     
     i_row( idx_sta : idx_end , 1 ) = dimC3 +         j_col_3;
     j_col( idx_sta : idx_end , 1 ) = dimC3 + n_123 + i_row_3;
     k_val( idx_sta : idx_end , 1 ) =       conj(k_val_3);

% - G_2' matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row(idx_sta : idx_end , 1 ) = dimC3 +             j_col_2;
     j_col(idx_sta : idx_end , 1 ) = dimC3 + 2 * n_123 + i_row_2;
     k_val(idx_sta : idx_end , 1 ) =             - conj(k_val_2);

% ======= Second block rows 

% - G_3' matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_3;
     
     i_row( idx_sta : idx_end , 1 ) = dimC3 + n_123 + j_col_3;
     j_col( idx_sta : idx_end , 1 ) = dimC3 +         i_row_3;
     k_val( idx_sta : idx_end , 1 ) =         - conj(k_val_3);
     
% - G_1 matrix         
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = dimC3 +     n_123 + i_row_1;
     j_col( idx_sta : idx_end , 1 ) = dimC3 + 2 * n_123 + j_col_1;
     k_val( idx_sta : idx_end , 1 ) =           - k_val_1;

% ===== Third block rows

% G_2' matrix
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_2;
     
     i_row( idx_sta : idx_end , 1 ) = dimC3 + 2 * n_123 + j_col_2;
     j_col( idx_sta : idx_end , 1 ) = dimC3 +             i_row_2;
     k_val( idx_sta : idx_end , 1 ) =           conj(k_val_2);
     
%   G_1 matrix     
     idx_sta                = idx_end +     1;   
     idx_end                = idx_end + len_1;
     
     i_row( idx_sta : idx_end , 1 ) = dimC3 + 2 * n_123 + i_row_1;
     j_col( idx_sta : idx_end , 1 ) = dimC3 +     n_123 + j_col_1;
     k_val( idx_sta : idx_end , 1 ) =             k_val_1;
%
% Final matrix C
%
     
     mtx_C            = sparse(i_row, j_col, k_val); 
     
     dimC4   = 4 * dimC; 
     
     i_row_A = [i_row; dimC4+j_col; 2*dimC4+j_col; 3*dimC4+i_row];
     j_col_A = [j_col; dimC4+i_row; 2*dimC4+i_row; 3*dimC4+j_col];
     k_val_A = [k_val; conj(k_val); conj(k_val); k_val];
     
     GEP_mtx_A = sparse(i_row_A, j_col_A, k_val_A);
 end
 
 %
 % ======================= C_1 matrix =================
 
 function [i_row,j_col,k_val] = Gamma_1(n_1,n_2,n_3,period_vec1,vec_k)
%  

% Gamma_1  Diagonal Block matrix G_1

%
tmp                =  2 * pi * 1i * (vec_k' * period_vec1);
mtx_G_1            = -eye( n_1 );
mtx_G_1( n_1 , 1 ) =  exp( tmp );
%
for ii = 2 : n_1
    mtx_G_1( ii - 1 , ii ) = 1;
end 

%  
[i_row_tmp,j_col_tmp,k_val_tmp] = find(mtx_G_1);
%
%  
n_23 = n_2 * n_3;
r_len                           = length(i_row_tmp);
i_row                           = zeros(r_len * n_23 , 1 );
j_col                           = zeros(r_len * n_23 , 1 );
k_val                           = zeros(r_len * n_23 , 1 );
for ii = 1:n_2*n_3
    jj                             = ( ii - 1 ) * n_1;
    kk                             = ( ii - 1 ) * r_len ;
    i_row( kk + 1 : kk + r_len ,1) = jj + i_row_tmp;
    j_col( kk + 1 : kk + r_len ,1) = jj + j_col_tmp;
    k_val( kk + 1 : kk + r_len ,1) =      k_val_tmp;
end

 end
 
 %
 % ======================= C_2 matrix =================
 
function [i_row,j_col,k_val] = Gamma_2(n_1,n_2,n_3,period_vec1,period_vec2,vec_k)
 
[i_row_tmp,j_col_tmp,k_val_tmp] = Generate_mtx_G_2(n_1,n_2,period_vec1,period_vec2,vec_k);

r_len                           = length( i_row_tmp );
i_row                           = zeros( r_len  * n_3 , 1 );
j_col                           = zeros( r_len  * n_3 , 1 );
k_val                           = zeros( r_len  * n_3 , 1 );

for ii = 1 : n_3
    jj                    = ( ii - 1 ) * n_1 * n_2;
    kk                    = ( ii - 1 ) *    r_len ;
    i_row( kk + 1 : kk + r_len , 1 ) = jj + i_row_tmp;
    j_col( kk + 1 : kk + r_len , 1 ) = jj + j_col_tmp;
    k_val( kk + 1 : kk + r_len , 1 ) =      k_val_tmp;
end

end

%
 function [i_row, j_col, k_val] = Generate_mtx_G_2(n_1,n_2,period_vec1,period_vec2,vec_k)
 
%  Gamma_2  Diagonal Block matrix G_2
     
% 
n_12 = n_1 * n_2;
i_row = zeros(2 * n_12 , 1);
j_col = zeros(2 * n_12 , 1);
k_val = zeros(2 * n_12 , 1);
% 
     
     for jj = 1 : n_2 - 1    
         
         ii                   = ( jj - 1 ) * n_1;
         
         kk                   = 2 * ii;
         
         i_row( kk + 1 : kk + n_1 , 1 ) = ii + 1 : ii + n_1;
         j_col( kk + 1 : kk + n_1 , 1 ) = ii + 1 : ii + n_1;
         k_val( kk + 1 : kk + n_1 , 1 ) =  -ones( n_1 , 1 );
         
         i_row( kk + n_1 + 1 : kk + 2 * n_1 , 1 ) = ii +       1 : ii +     n_1;
         j_col( kk + n_1 + 1 : kk + 2 * n_1 , 1 ) = ii + n_1 + 1 : ii + 2 * n_1;
         k_val( kk + n_1 + 1 : kk + 2 * n_1 , 1 ) =             ones( n_1 , 1 );
     
     end     
     
 %  
     
     ii                   = ( n_2 - 1 ) * n_1;
     kk                   =           2 *  ii;
     
     i_row( kk + 1 : kk + n_1 , 1 ) = ii + 1 : ii + n_1;
     j_col( kk + 1 : kk + n_1 , 1 ) = ii + 1 : ii + n_1;
     k_val( kk + 1 : kk + n_1 , 1 ) =  -ones( n_1 , 1 );
     
     %! The inner product in Diamond case does not equal to length_y * vec_k(2) * 1i
     %! We should compute this term by this formula          
     %!                                          KdP = period_vec2 * vec_k;
     %!                                          tmp =         KdP *    li; 
     %! End
     
     KdP_1 = vec_k' * period_vec1;
     KdP_2 = vec_k' * period_vec2;     
     tmp_1 = 2 * pi * KdP_1 * 1i;
     tmp_2 = 2 * pi * KdP_2 * 1i;
    
     %  
     i_row( kk + n_1 + 1 : kk + (3/2) * n_1 , 1 ) =          ii + 1 : ii + (1/2) * n_1;
     j_col( kk + n_1 + 1 : kk + (3/2) * n_1 , 1 ) = (1/2) * n_1 + 1 :              n_1;
     k_val( kk + n_1 + 1 : kk + (3/2) * n_1 , 1 ) = exp( tmp_2 - tmp_1 ) * ones( (1/2) * n_1 , 1 );
     
     %  
     i_row( kk + (3/2) * n_1 + 1 : kk + 2 * n_1 , 1 ) = ii + (1/2) * n_1 + 1 :    ii + n_1;
     j_col( kk + (3/2) * n_1 + 1 : kk + 2 * n_1 , 1 ) =                    1 : (1/2) * n_1;
     k_val( kk + (3/2) * n_1 + 1 : kk + 2 * n_1 , 1 ) = exp(tmp_2)         * ones( (1/2) * n_1 , 1 );
 end
 
 %
 % ======================= C_3 matrix =================
  function [i_row, j_col, k_val] = Gamma_3(n_1,n_2,n_3,period_vec1,period_vec2,period_vec3,vec_k)
 
%         
     n_12   = n_1  * n_2;
     n_123  = n_12 * n_3;
     i_row  = zeros( 2 * n_123 , 1 );
     j_col  = zeros( 2 * n_123 , 1 );
     k_val  = zeros( 2 * n_123 , 1 );
     
 %     
     for jj = 1 : n_3 - 1
         
         ii                   = ( jj - 1 ) * n_12;
         
         kk                   =          2 *   ii;
         
         i_row( kk + 1 : kk + n_12 , 1 ) = ii + 1 : ii + n_12;
         j_col( kk + 1 : kk + n_12 , 1 ) = ii + 1 : ii + n_12;
         k_val( kk + 1 : kk + n_12 , 1 ) =  -ones( n_12 , 1 );
                  
         i_row( kk + n_12 + 1 : kk + 2 * n_12 , 1 ) = ii + 1 : ii + n_12;
         j_col( kk + n_12 + 1 : kk + 2 * n_12 , 1 ) =          ii + n_12 + 1 : ii + 2 * n_12;
         k_val( kk + n_12 + 1 : kk + 2 * n_12 , 1 ) =                       ones( n_12 , 1 );
     end
  
% ( Simple Cube )  
%     ii                   = ( n_3 - 1 ) * n_12;
%     kk                   =           2 *   ii;
%
%     i_row( kk + 1 : kk + n_12 , 1 ) = ii + 1 : ii + n_12;
%     j_col( kk + 1 : kk + n_12 , 1 ) = ii + 1 : ii + n_12;
%     k_val( kk + 1 : kk + n_12 , 1 ) =  -ones( n_12 , 1 );
   
%     tmp                                       = length_z * vec_k(3) * 1i;    
%    i_row( kk + n_12 + 1 : kk + 2 * n_12 , 1 ) = ii + 1 : ii + n_12;
%    j_col( kk + n_12 + 1 : kk + 2 * n_12 , 1 ) =      1 :      n_12;
%    k_val( kk + n_12 + 1 : kk + 2 * n_12 , 1 ) =  exp( tmp ) * ones( n_12 , 1 ); 
  
% ( Diamond Structrue )

      ii                   = ( n_3 - 1 ) * n_12;
      kk                   =           2 *   ii;

      i_row( kk + 1 : kk + n_12 , 1 ) = ii + 1 : ii + n_12;
      j_col( kk + 1 : kk + n_12 , 1 ) = ii + 1 : ii + n_12;
      k_val( kk + 1 : kk + n_12 , 1 ) =  -ones( n_12 , 1 );
     
     %! The inner product in diamond case does not equal to length_z * vec_k(3) * 1i
     %! We should compute this term by this formula          
     %!                                                      KdP = period_vec3 * vec_k;
     %!                                                      tmp =            KdP * 1i; 
     %! End 
     
      KdP_2 = vec_k' * period_vec2;
      KdP_3 = vec_k' * period_vec3;                                                    
      tmp_2 = 2 * pi * KdP_2 * 1i;
      tmp_3 = 2 * pi * KdP_3 * 1i;
      
      %  
      i_row( kk + n_12 + 1 : kk + (4/3) * n_12 , 1 ) = ii + 1 : ii + (1/3) * n_12;      
      j_col( kk + n_12 + 1 : kk + (4/3) * n_12 , 1 ) =               (2/3) * n_12 + 1 : n_12;
      k_val( kk + n_12 + 1 : kk + (4/3) * n_12 , 1 ) = exp( tmp_3 - tmp_2 ) * ones( (1/3) * n_12 , 1 ); 
      
      %  
      
      % ?c?y (2/3) * n_12 x n_12 block G_2 ?x?}       
      
      [row_tmp, col_tmp, val_tmp] = Generate_mtx_G_2(n_1,n_2,period_vec1,period_vec2,vec_k);
     
      val_tail_1 = exp(tmp_3-tmp_2) * val_tmp( 2 * n_12 - n_1 + 1 : 2 * n_12 - (1/2) * n_1 );                                   
      val_tail_2 = exp(tmp_3-tmp_2) * val_tmp(                      2 * n_12 - (1/2) * n_1 + 1 : 2 * n_12 );        
      
     
      idx_start = kk + (4/3) * n_12;
      idx_end   = kk + (4/3) * n_12 + (1/2) * n_1;
      
      idx_row   = n_123 - (2/3) * (n_12);
      idx_col   =         (1/2) * n_1;
      
      for iii = 1 : (2/3) * n_2 
          
      i_row( idx_start + 1 : idx_end, 1 ) = idx_row + 1 : idx_row + (1/2) * n_1;
      j_col( idx_start + 1 : idx_end, 1 ) = idx_col + 1 : idx_col + (1/2) * n_1;
      k_val( idx_start + 1 : idx_end, 1 ) =                          val_tail_1;
      
      idx_start = idx_end;
      idx_end   = idx_end + (1/2) * n_1;
      idx_row   = idx_row + (1/2) * n_1;
      idx_col   = idx_col - (1/2) * n_1;
      
      i_row( idx_start + 1 : idx_end, 1 ) = idx_row + 1 : idx_row + (1/2) * n_1;
      j_col( idx_start + 1 : idx_end, 1 ) = idx_col + 1 : idx_col + (1/2) * n_1;
      k_val( idx_start + 1 : idx_end, 1 ) =                          val_tail_2;
      
      idx_start = idx_end;
      idx_end   = idx_end + (1/2) * n_1;
      idx_row   = idx_row + (1/2) * n_1;
      idx_col   = idx_col + (3/2) * n_1;
      
      end
     
  end
