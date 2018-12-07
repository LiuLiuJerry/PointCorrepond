function [G,Ceq]=nonlcon_imp(X, xyz_kpts, nbones)
%X : n *(nbones * 3) 2维数组
%点到线段的距离最短
n = size(X, 1);
% nbones = 4;
len_X = zeros(n,nbones-1);
len_pts = zeros(n,1);
v_last = zeros(n, 3);
vl_x = zeros(n, 3);
% X = reshape(X, [size(xyz_pts,1), nbones*3]);
for i = 1:n
    xyz = X(i,:);
    pts = xyz_kpts{i};
    % 骨节长度和约束 G1
    diff_pts = diff(pts);
    dis_p2p = sqrt(dot(diff_pts, diff_pts, 2));
    len_pts(i,:) = sum(dis_p2p);
    % 距离等长约束 Seq, 不考虑固定点
    xyz = reshape(xyz, [nbones,3]);
    diff_X = diff(xyz);
    len_X(i, :) = sqrt(dot(diff_X, diff_X, 2))'; %各个骨头的长度
    % G2
    v_last(i,:) = pts(end,:) - pts(end-1,:);
    vl_x(i,:) = xyz(end,:) - pts(end-1,:);
end
% 估算的骨节的长度之和应该比采样点的距离之和小
% G = sum(len_X,2) - max(len_pts);
G = -1.*dot(v_last, vl_x, 2);
% 各个腿对应骨节长度相等
len_X = [len_X(end,:); len_X];
diff_len = abs(diff(len_X));
Ceq = sum(diff_len, 2);
% Ceq = sum(diff_len)';
