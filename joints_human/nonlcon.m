function [G,Ceq]=nonlcon(X, xyz_pts, nbones)
%X : n *(nbones * 3) 2ά����
%�㵽�߶εľ������
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
    % �ǽڳ��Ⱥ�Լ�� G1
    diff_pts = diff(pts);
    dis_p2p = sqrt(dot(diff_pts, diff_pts, 2));
    len_pts(i,:) = sum(dis_p2p);
    % ����ȳ�Լ�� Seq
    startpts = pts(1,:);
    xyz = reshape(xyz, [nbones,3]);
    xyz = [startpts; xyz];
    diff_X = diff(xyz);
    len_X(i, :) = sqrt(dot(diff_X, diff_X, 2))'; %������ͷ�ĳ���
    % G2
    v_last(i,:) = pts(end,:) - pts(end-1,:);
    vl_x(i,:) = xyz(end,:) - pts(end-1,:);
end
% ����Ĺǽڵĳ���֮��Ӧ�ñȲ�����ľ���֮��С
G1 = sum(len_X,2) - max(len_pts);
% G2 = [];
G2 = -1.*dot(v_last, vl_x, 2); %���һ���ؼ��������һ�������㾡���ܽ�
G = [G1;G2];
% �����ȶ�Ӧ�ǽڳ������
len_LR = abs(len_X(1:end/2, :) - len_X(end/2+1:end, :));
diff_len = zeros(n, n);
for i = 1:n
    % ģ��֮���������
    dis = abs(len_X(i,:) - len_X);% n*nbones
    diff_len(i,:) = sum(dis, 2)';
end
diff_len = sum(diff_len,2);
Ceq1 = sum(len_LR, 2);
Ceq2 = sum(diff_len, 2); 
Ceq = [Ceq1;Ceq2];
