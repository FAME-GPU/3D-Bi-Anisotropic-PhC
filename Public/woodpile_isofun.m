function Point_idx_tmp = woodpile_isofun(Point_set,w,a)
% 判断点集Point_set是否处于介质内
% w为长方体宽度,a为晶格常数
x = Point_set(:,1);
y = Point_set(:,2);
z = Point_set(:,3);
t1 = y-x;  t2 = y+x;  w = sqrt(2)/2*w;
Point_idx_tmp = (   (z>=0 & z<a/4)   & ( (t1>=(a-w)) | (t1<=-(a-w)) | ((t1>=-w) & (t1<=w))  )  )...
              | ( (z>=a/4 & z<a/2)   & ( ((t2>=3*a/2-w) & (t2<=3*a/2+w)) | ((t2>=a/2-w) & (t2<=a/2+w))  )  )...
              | ( (z>=a/2 & z<3*a/4) & ( ((t1>=(-a/2-w))  & (t1<=(-a/2+w)))  | ((t1>=(a/2-w)) & (t1<=(a/2+w)))  )  )...
              | ( (z>=3*a/4 & z<=a)  & ( (t2>=(2*a-w)) | (t2<=w) | ((t2>=a-w) & (t2<=a+w))  )  );
Point_idx_tmp = double(Point_idx_tmp);
Point_idx_tmp = find(Point_idx_tmp~=0);

% data = Point_set([1;Point_idx;216000],:);
% x1 = data(:, 1);
% y1 = data(:, 2);
% z1 = data(:, 3);
% figure;
% scatter3(x1, y1, z1, 'filled'); % 绘制填充的散点图
% title('3D Scatter Plot of 100 Points');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% grid on;
% scatter3(x1, y1, z1, 36, z1, 'filled');  % 36 是点的大小，z 是颜色数据
% colorbar;  % 显示颜色条

end