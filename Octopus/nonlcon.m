function [G,Ceq]=nonlcon(X, xyz_pts, nbones)
%X : n *(nbones * 3) 2维数组
%点到线段的距离最短
n = size(X, 1);
% nbones = 4;
len_X = zeros(n,nbones);
len_pts = zeros(n,1);
v_last = zeros(n, 3);
vl_x = zeros(n, 3);
% X = reshape(X, [size(xyz_pts,1), nbones*3]);
for i = 1:n
    xyz = X(i,:);
    pts = xyz_pts{i};
    % 骨节长度和约束 G1
    diff_pts = diff(pts);
    dis_p2p = sqrt(dot(diff_pts, diff_pts, 2));
    len_pts(i,:) = sum(dis_p2p);
    % 距离等长约束 Seq
    startpts = pts(1,:);
    xyz = reshape(xyz, [nbones,3]);
    xyz = [startpts; xyz];
    diff_X = diff(xyz);
    len_X(i, :) = sqrt(dot(diff_X, diff_X, 2))'; %各个骨头的长度
    % G2
    v_last(i,:) = pts(end,:) - pts(end-1,:);
    vl_x(i,:) = xyz(end,:) - pts(end-1,:);
end
% 估算的骨节的长度之和应该比采样点的距离之和小
G1 = sum(len_X,2) - max(len_pts);
% G2 = [];
G2 = -1.*dot(v_last, vl_x, 2); %最后一个关键点离最后一个采样点尽可能近
G = [G1;G2];
% 各个腿对应骨节长度相等
len_LR = abs(len_X(1:end/2, :) - len_X(end/2+1:end, :));
diff_len = zeros(n, n);
for i = 1:n
    % 模型之间两两相等
    dis = abs(len_X(i,:) - len_X);% n*nbones
    diff_len(i,:) = sum(dis, 2)';
end
diff_len = sum(diff_len,2);
Ceq1 = sum(len_LR, 2);
Ceq2 = sum(diff_len, 2); 
Ceq = [Ceq1;Ceq2];
