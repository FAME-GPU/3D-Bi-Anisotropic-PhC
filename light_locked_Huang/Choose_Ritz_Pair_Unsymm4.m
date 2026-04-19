function [ewno, ritz_ew, ritz_ev] = Choose_Ritz_Pair_Unsymm4(iteno, ew, ew_order, target, VRR, ...
    vec_r_tol, Choose_Type, RestartProjProbDim, target_org, Choose_Type_org)

% --- select the desired eigenpair (lambda,x)
%     such that lambda is the
%     closest to target

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   New version, setup at Aug. 24, 2004
%
global zflag_ew flag_no_ew
global flag_gpu

warning off all

sort_ew     = zeros(iteno,1);
ewp_idx     = zeros(iteno,1);
sort_No     = 5;
if flag_gpu
    sort_ew = gpuArray(sort_ew);
    ewp_idx = gpuArray(ewp_idx);
end

sort_ritz = ew;

if ( ew_order > 1 )
    
    mult_cn_ew    = zeros(flag_no_ew,1);
    mult_ritz_w   = zeros(iteno,1);
    ld            = zeros(iteno,1);
    idx_mult_ritz = zeros(iteno,1);
    if flag_gpu
        mult_cn_ew    = gpuArray(mult_cn_ew);
        mult_ritz_w   = gpuArray(mult_ritz_w);
        ld            = gpuArray(ld);
        idx_mult_ritz = gpuArray(idx_mult_ritz);
    end
    
    sort_ew(1:flag_no_ew) = zflag_ew(1:flag_no_ew);
    
    %
    % mult_cn_ew(i) : multiplicity of eigenvalue flag_ew(i)
    %
    % Find the multiplicity of eigenvalue. Suppose that
    %        \lambda_1 = ... = \lambda_k
    % are mutiple eigenvalues. Then
    %        mult_cn_ew(1) = k
    % and
    %        mult_cn_ew(j) = 0, for j = 2, ..., k
    %
    for i = 1:flag_no_ew   % 判断已经计算出的特征值sort_ew的重根数mult_cn_ew     flag_no_ew为已经计算出的特征值个数，iteno为ritz值个数
        if ( real(sort_ew(i)) > - 1.0e13 )
            no_cn_ew  = 1;
            for i1 = i + 1:flag_no_ew
                if ( real(sort_ew(i1)) > - 1.0e13 )
                    k        = -1;
                    flag_tol =  1;
                    
                    while ( k <= 1 && flag_tol == 1 )
                        if ( abs(sort_ew(i1)-sort_ew(i)) <= 10.0^(k)*vec_r_tol  ) % 判断规则：sort_ew(i1)-sort_ew(i) < =0.1 or 1 or 10 * 1e-9
                            flag_tol = 0;
                        else
                            k = k + 1;
                        end
                    end
                    
                    if ( flag_tol == 0 )
                        no_cn_ew    = no_cn_ew + 1;
                        sort_ew(i1) = -1.0e14;  % 将第二个及以后的重根置为-1.0e14，则其不满足 > - 1.0e13，不再参与后续判断
                    end
                end
            end
        else
            no_cn_ew = 0;
        end
        mult_cn_ew(i) = no_cn_ew;
    end
    
    %
    % mult_ritz_w(i) : multiplicity of ritz value ALPHAR(i)
    % idx_ritz_w(j)  : denote the ritz value ALPHAR(j) to
    %                  be equal to ritz value ALPHAR(idx_ritz_w(j))
    %
    % Find the multiplicity of ritz value. Suppose that
    %        \theta_1 = ... = \theta_k
    % are mutiple ritz values. Then
    %        mult_ritz_w(1) = k
    % and
    %        mult_ritz_w(j) = 0, for j = 2, ..., k
    %
    
    idx_ritz_w = zeros(iteno,1);
    if flag_gpu
        idx_ritz_w = gpuArray(idx_ritz_w);
    end
    for i = 1:iteno   % 1.判断已经计算出的ritz值的重根数no_ritz_w 与上面相同的判断规则 2.给出ritz值的索引idx_ritz_w 规则为首次出现的ritz值索引0，其他次出现的ritz值索引为首次出现的相同ritz值的位置号
        if ( idx_ritz_w(i) == 0 )
            no_ritz_w = 1;
            for i1 = i + 1:iteno
                if ( idx_ritz_w(i1) == 0 )
                    k        = -1;
                    flag_tol =  1;
                    while ( k <= 1 && flag_tol == 1 )
                        if ( abs(sort_ritz(i1)-sort_ritz(i)) <= 10.0^(k)*vec_r_tol  )
                            flag_tol = 0;
                        else
                            k = k + 1;
                        end
                    end
                    if ( flag_tol == 0 )
                        no_ritz_w      = no_ritz_w + 1;
                        idx_ritz_w(i1) = i; % 比如有ritz值[2 2 3 2 3] 则算出的idx_ritz_w = [0 1 0 1 3]
                    end
                end
            end
        else
            no_ritz_w = 0;
        end
        mult_ritz_w(i) = no_ritz_w;
    end
    %
    % The ritz values which are equal to the convergent eigenvalues
    % are set to be infty.
    %
    for i = 1:flag_no_ew   % 判断ritz是否与已收敛的特征值相等，对每个特征值，都遍历所有ritz值
        if ( real(sort_ew(i)) > - 1.0e13 )
            i1   = 1;
            flag = 1;
            while ( i1 <= iteno && flag == 1 )
                if ( idx_ritz_w(i1) == 0 )   % 对多次出现的相同ritz值，只需观察第一次出现的ritz值即可
                    k        = -1;
                    flag_tol = 1;
                    while ( k <= sort_No && flag_tol == 1 )
                        if ( abs(sort_ritz(i1)-sort_ew(i)) <= 10.0^(k) * vec_r_tol  )% 判断规则：sort_ritz(i1)-sort_ew(i) < =0.1 1 10 1e2 1e3 1e4 1e5 * 1e-9
                            flag_tol = 0;
                            ind_sort = i1;  % ind_sort：与当前特征值相同的ritz值的索引
                        else
                            k = k + 1;
                        end
                    end
                    if ( flag_tol == 0 )
                        flag = 0;
                    else
                        i1   = i1 + 1;
                    end
                else
                    i1 = i1 + 1;
                end
            end
            
            if ( flag == 1 ) % 若现有没有能与已收敛的特征值相等的ritz值，报以下错：
                fprintf('Error in Choose_Ritz_Pair\n');
                % fprintf('cn_ew(%3.0f) = \n', i, zflag_ew(i));
                fprintf('cn_ew(%d) = %f+%fi\n', i, real(zflag_ew(i)), imag(zflag_ew(i)));
                fprintf('Convergent eigenvalues are:\n');
                disp(zflag_ew)
                fprintf('\n');
                fprintf('Ritz values are:\n');
                disp(sort_ritz)
                fprintf('\n');
            else
                if ( mult_cn_ew(i) == mult_ritz_w(ind_sort) ) % 若第i个特征值的重数等于第ind_sort个ritz值的重数
                    sort_ritz(ind_sort) = 1.0e14;
                    for jj = 1:iteno
                        if ( idx_ritz_w(jj) == ind_sort )
                            sort_ritz(jj) = 1.0e14;
                        end
                    end
                    idx_ritz_w(ind_sort) = -1;
                elseif ( mult_cn_ew(i) < mult_ritz_w(ind_sort) )  % 若第i个特征值的重数 < 第ind_sort个ritz值的重数
                    %
                    % Another multiple eigenvalue which does not contain
                    % in the set of convergent eigenvalues belongs to
                    % the set of ritz values.
                    % The ritz values which are equal to the convergent
                    % multiple eigenvalues are set to be infty so that
                    % the corresponding eigenvector does not be chosen.
                    % How to choose these ritz values? Let span(V1) be the
                    % eigenspace associated with the convergent eigenvalues
                    % and V = [ V1 V2 ], s = [ s1^T  s2^T ]^T. Then
                    %               u = V s
                    % is an approximated eigenvector. If s2 = 0, then u
                    % is belong to span(V1) which is an already convergent
                    % eigenvector. Hence u is not a target vector. We want
                    % to choose u such that u is not belong to span(V1) and
                    % {V1, u} is linearly independent.
                    %
                    k                = 1;
                    idx_mult_ritz(k) = ind_sort;                      % 记录第一个重复里兹值的索引
                    ld(k)            = norm(VRR(ew_order:iteno,ind_sort));
                    for jj = 1:iteno
                        if ( idx_ritz_w(jj) == ind_sort )             % 找到所有与ind_sort重复的里兹值
                            k                = k + 1; 
                            idx_mult_ritz(k) = jj;                    % 记录所有重复里兹值的索引
                            ld(k)            = norm( VRR(ew_order:iteno,jj) );  % 计算每个重复里兹向量的模长
                        end
                    end
                    if ( abs(k - mult_ritz_w(ind_sort)) > 0 )
                        fprintf('error in multiplicity of ritz value \n');
                    end
                    for jj = 1:mult_cn_ew(i)
                        tmp = ld(1);
                        k   = 1;
                        for i1 = 2:mult_ritz_w(ind_sort)
                            if ( tmp > ld(i1) )
                                tmp = ld(i1);                         % 找到模长最小的重复的里兹向量
                                k   = i1;                             % 模长最小的重复的里兹向量对应的索引
                            end
                        end
                        sort_ritz(idx_mult_ritz(k)) = 1.0e14;
                        ld(k)                       = 1.0e14;
                    end
                    idx_ritz_w(ind_sort) = -1;
                else                                      % 若第i个特征值的重数 > 第ind_sort个ritz值的重数
                    sort_ritz(ind_sort) = 1.0e14;         % 标记所有相关里兹值为无效
                    for jj = 1:iteno
                        if ( idx_ritz_w(jj) == ind_sort )
                            sort_ritz(jj) = 1.0e14;
                        end
                    end
                    idx_ritz_w(ind_sort) = -1;
                    fprintf('Error for multiplicity of convergent \n');
                    fprintf('eigenvalue > multiplicity of ritz value \n');
                end
            end
        end
    end
end

%
% Choose the ritz values which are closest the target value or
% approximate to the convergent eigenvalues.
%
jj          = 1;              % 计数器：记录符合条件的里兹值数量+1
dist_target = zeros(iteno,1); % 存储里兹值与目标值的距离
if flag_gpu
    dist_target = gpuArray(dist_target);
end

switch Choose_Type
    case('CL')
        %      choose the ritz value which is the closest target value
        %
        for i = 1:iteno
%             if ( abs(sort_ritz(i)) < 1.0e10 )
%                 ind_sort = 0;
%                 for k = 1:flag_no_ew
%                     if ( abs(sort_ritz(i)-zflag_ew(k)) < 1.0e3*vec_r_tol)
% %                     if ( abs(abs(sort_ritz(i))-abs(zflag_ew(k))) < 1.0e3*vec_r_tol) % ss-coorect
%                         ind_sort = k;
%                     end
%                 end
%                 if ( ind_sort > 0 )
%                     dist_target(jj) = abs(sort_ritz(i)-zflag_ew(ind_sort));
% %                     dist_target(jj) = abs( abs(sort_ritz(i))-abs(zflag_ew(ind_sort)) ); % ss-correct
%                 else
%                     dist_target(jj) = abs(sort_ritz(i)-target);
% %                     dist_target(jj) = abs( abs(sort_ritz(i))-abs(target) ); % ss-correct
%                 end
%                 ewp_idx(jj) = i;
%                 jj          = jj + 1;
%             end

            % ss-correct
            if ( abs(sort_ritz(i)) < 1.0e10 )
                %if real(ew(i)) > real(target)
                if ( abs(real(ew(i))) < 1e-4 && imag(ew(i)) > 1e-3 ) || ( real(ew(i)) > real(target) && abs(imag(ew(i))) < 1e-8 )
                %if ( real(ew(i)) > real(target) || abs(imag(ew(i))) > 1e-6 ) %--- && abs(imag(ew(i))) < 1e-6 )---|| abs(imag(ew(i))) > 1e-6 )
                %if ( abs(real(ew(i))) < 0.2 && imag(ew(i)) > 1e-6 )
                    ewp_idx(jj)     = i;       % ewp_idx: 满足条件的ritz值的索引
                    dist_target(jj) = abs(real(ew(i))-real(target));
                    jj              = jj + 1;  % jj: 满足条件的ritz值的个数+1
                end
            end
        end
    case('RGTR')
        %      if the target value is 'Real' and the wanted eigenvalue is
        %      'Greater' 'Than' target. The Ritz value in every iteration will
        %      be 'Real'
        for i = 1:iteno
            if ( abs(sort_ritz(i)) < 1.0e10 )
                if ( real(ew(i)) > real(target) && abs(imag(ew(i))) < 1.0e-12 )
                    ewp_idx(jj)     = i;
                    dist_target(jj) = abs(real(ew(i))-real(target));
                    jj              = jj + 1;
                end
            end
        end
    case('RGTC')
        %      if the target value is 'Real' and the real part of the wanted
        %      eigenvalue is 'Greater' 'Than' target. The Ritz value in every
        %      iteration can be 'Complex'
        for i = 1:iteno
            if ( abs(sort_ritz(i)) < 1.0e10 )
                if ( real(ew(i)) > real(target) )
                    ewp_idx(jj)     = i;
                    dist_target(jj) = abs(ew(i)-target);
%                     dist_target(jj) = abs(abs(ew(i))-target); % ss-correct
                    jj              = jj + 1;
                end
            end
        end
    otherwise
        fprintf('In-correct Choose_Type; retype again \n');
end

if ( jj == 1 )      % 若没有满足条件的ritz值
    for i = 1:iteno
        ewp_idx(i)     = i;   % ewp_idx: 满足条件的ritz值的索引，此行为取所有索引
        dist_target(i) = abs(sort_ritz(i)-target);    % dist_target：满足条件的ritz值的实部与目标值的距离
%         dist_target(i) = abs( abs(sort_ritz(i))-target ); % ss-coorect
    end
    jj = iteno + 1;
end

ewno = min(jj-1, RestartProjProbDim);   % 选择ritz对的数量不超过重启投影空间的维数

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   New version, setup at Jan. 09, 2004
%
ritz_ev = zeros(iteno,ewno);
ritz_ew = zeros(ewno,1);
idx     = zeros(RestartProjProbDim,1);
if flag_gpu
    ritz_ev = gpuArray(ritz_ev);
    ritz_ew = gpuArray(ritz_ew);
    idx = gpuArray(idx);
end

for i = 1:ewno
    ind_sort = 1;
    for k = 2:jj-1              % 选出与目标值距离最小的ritz值的索引
        if ( dist_target(k) < dist_target(ind_sort) )
            if ( abs(target_org - target) < 1.0e-10 )
%             if ( abs(abs(target_org) - target) ) % ss-correct
                ind_sort = k;
            else
                switch Choose_Type_org
                    case 'CL'
                        ind_sort = k;
                    case 'RGTR'
                        if ( real(ew(ewp_idx(k))) > real(target_org) && ...
                                abs(imag(ew(ewp_idx(k)))) < 1.0e-12 )
                            ind_sort = k;
                        end
                    case 'RGTC'
                        if ( real(ew(ewp_idx(k))) > real(target_org) )
                            ind_sort = k;
                        end
                end
            end
        end
    end
    
    if ( dist_target(ind_sort) > 1.0e12 )
        [~,ind_sort] = min( dist_target(1:jj-1) );
    end
    
    idx(i,1)              = ewp_idx(ind_sort);
    ritz_ew(i)            = ew(ewp_idx(ind_sort));
    ritz_ev(1:iteno,i)    = VRR(1:iteno,ewp_idx(ind_sort));
    ritz_ev(1:iteno,i)    = ritz_ev(1:iteno,i) / norm(ritz_ev(1:iteno,i));
    dist_target(ind_sort) = 1.0e14;
end

%
%   modify at 2013/08/21
%
if ( ewno < RestartProjProbDim )
    out_idx = zeros(iteno-ewno,1);          % out_idx：未被选中的ritz对索引
    if flag_gpu
        out_idx = gpuArray(out_idx);
    end
    kk      = 1;
    for ii = 1:iteno
        if ( isempty(find( abs(ii-idx) == 0, 1 )) )
            out_idx(kk) = ii;
            kk          = kk + 1;
        end
    end
    
    %      if the target value is 'Real' and the wanted eigenvalue is
    %      'Greater' 'Than' target. The Ritz value in every iteration will
    %      be 'Real'
    jj = 1;                                 % 从剩余里兹对中筛选候选（更新 ewp_idx 和 dist_target）
    for i = 1:iteno 
        if ( abs(sort_ritz(i)) < 1.0e10 )
            if ( ~isempty(find( abs(i-out_idx) == 0, 1 )) )
                %if ( real(ew(i)) > real(target) )
                ewp_idx(jj)     = i;
                dist_target(jj) = abs(ew(i) - target);           % 计算该里兹值与目标值的距离用于排序
%                 dist_target(jj) = abs(abs(ew(i)) - target); % ss-correct
                jj              = jj + 1;                        % 更新ritz对数量
                %end
            end
        end
    end
    
    [~,idx2] = sort(dist_target(1:jj-1));                         % 对补充候选的距离进行排序，得到排序索引idx2（距离越小越优先）
    if ( ~isempty(idx2) )
        for ii = 1:min(RestartProjProbDim-ewno, length(idx2))     % 按排序结果依次补充里兹值和里兹向量
            ritz_ew(ewno+ii)         = ew(ewp_idx(idx2(ii)));
            ritz_ev(1:iteno,ewno+ii) = VRR(1:iteno,ewp_idx(idx2(ii)));
            ritz_ev(1:iteno,ewno+ii) = ritz_ev(1:iteno,ewno+ii) / ...
                norm(ritz_ev(1:iteno,ewno+ii));
        end
        ewno = ewno + min(RestartProjProbDim-ewno, length(idx2));
    end
end

