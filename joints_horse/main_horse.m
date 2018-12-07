clear all;
close all;
%% 骨架和提取骨架后得到的矩阵啊之类的
formatSpec = 'D:/GitHub/data/Chen/Skeletons/%s/mesh_%.4d.ply';
P_formatSpec = 'D:/GitHub/data/Chen/Skeletons/%s/mesh_%.4d_skeleton';
idx = [0:47];

% 定义骨架的结构
n_imptpts = 3;
nbones_h = 1;
nbones_leg = 3;
nbones_t = 3; 

% 声明数据
N = size(idx, 2);
imptpts = cell(N, 1);
ear = cell(N,1);
head = cell(N, 1);
legsf_r = cell(N, 1);
legsb_r = cell(N, 1);
legsf_l = cell(N, 1);
legsb_l = cell(N, 1);
tail = cell(N, 1);
backbone = cell(N,1);

cnt_models = 0;
xyz = cell(N, 1);
xyz_models = cell(N,1);
for ii = idx
paths = sprintf(formatSpec, 'horse_off',ii);
P_path = sprintf(P_formatSpec, 'horse_off', ii);
ptCloud = pcread(paths);
P = load(P_path);
P = P.P;
bbox = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
rs = bbox(:,1)-bbox(:,2);
diameter = sqrt(dot(rs,rs));

% init key points
[importantpts, parts, legslength]  = keypts_horse(ptCloud.Location, P, 0);
cnt_models = cnt_models+1;

xyz{cnt_models} = ptCloud.Location;
xyz_model = ptCloud.Location;
xyz_models{cnt_models} = xyz_model;
% 直接存储坐标
imptpts{cnt_models} = xyz_model(importantpts, :);%z轴从小到大排序
% ear{cnt_models} = xyz_model(fliplr(parts{1}), :);
head{cnt_models} = xyz_model(fliplr(parts{1}), :);
% (-1,-1),(1,-1),(-1,1),(1,1)右后，左后，右前，左前
legsb_r{cnt_models} = xyz_model(fliplr(parts{2}), :); %反转向量，从里往外
legsf_r{cnt_models} = xyz_model(fliplr(parts{3}), :);
legsb_l{cnt_models} = xyz_model(fliplr(parts{4}), :);
legsf_l{cnt_models} = xyz_model(fliplr(parts{5}), :);
tail{cnt_models} = xyz_model(fliplr(parts{6}), :);
backbone{cnt_models} = xyz_model(parts{7}, :);
end

%% 优化关键点
% 用开始提取的关键点做x0
x0_kpts = zeros(cnt_models, n_imptpts*3);
for i = 1:cnt_models
    x0_tmp = imptpts{i}([1,floor(end/2)+1,end]', :);
    x0_kpts(i, :) = x0_tmp(:);
end
% 先优化脊椎骨的关节点
Xkpts = fmincon(@(X)fun_impdist(X, backbone, imptpts, n_imptpts), x0_kpts, [], [], [], [], [], [], @(X)nonlcon_imp(X,  imptpts, n_imptpts));

%% 优化骨节连接点
% 估算x0
xh0 = zeros(cnt_models, nbones_h*3);
xf0 = zeros(cnt_models*2, nbones_leg*3);
xb0 = zeros(cnt_models*2, nbones_leg*3);
xt0 = zeros(cnt_models, nbones_t*3);
for t = 1:cnt_models
    kpts = reshape(Xkpts(t, :), [n_imptpts, 3]);
    kpts = sortrows(kpts,3);
    head{t} = [kpts(3,:); head{t}];  
    legsf_l{t} = [kpts(2,:);legsf_l{t}];
    legsf_r{t} = [kpts(2,:);legsf_r{t}];
    legsb_l{t} = [kpts(1,:);legsb_l{t}];
    legsb_r{t} = [kpts(1,:);legsb_r{t}];    
    tail{t} = [kpts(1,:);tail{t}];
end
    
for t = 1:cnt_models
    init_h = head{t}(end, :);
    xh0(t, :) = init_h(:);
    initf_r = legsf_r{t}([floor(end/2), end-1:end]', :);
    initf_l = legsf_l{t}([floor(end/2), end-1:end]', :);
    xf0(t,:) = initf_l(:); %先排n个左腿，再排n个右腿
    xf0(cnt_models+t, :) = initf_r(:);
    initb_r = legsb_r{t}([floor(end/2), end-1:end]', :);
    initb_l = legsb_l{t}([floor(end/2), end-1:end]', :);
    xb0(t,:) = initb_l(:);
    xb0(cnt_models+t, :) = initb_r(:);
    size_t = floor(size(tail{t}, 1)/3);
    init_t = tail{t}([size_t, size_t*2, end]', :);
    xt0(t,:) = init_t(:);
end

VLB = [];
VUB = [];
% 应该把关节点也加入坐标中
legsf = [legsf_l, legsf_r];
legsb = [legsb_l, legsb_r];
% 传入坐标以前腿的坐标应该减去最开始的坐标
% 优化腿和头的关节点
Xh = fmincon(@(X)fun_dist(X, head, nbones_h), xh0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, head, nbones_h));
Xf = fmincon(@(X)fun_dist(X, legsf, nbones_leg), xf0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, legsf, nbones_leg));
Xb = fmincon(@(X)fun_dist(X, legsb, nbones_leg), xb0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, legsb, nbones_leg));
Xt = fmincon(@(X)fun_dist(X, tail, nbones_t), xt0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, tail, nbones_t));
%% plot
for ii = 1:cnt_models
    figure(idx(ii)+100);
    xyz_model = xyz_models{ii}; 
    scatter3(xyz_model(:,1), xyz_model(:,2), xyz_model(:,3), 15, 'filled'); hold on;

    % head
    xyz_heads = reshape(Xh(ii,:), [nbones_h, 3]);
    xyz_heads = [head{ii}(1,:); xyz_heads];
    plot3(xyz_heads(:, 1), xyz_heads(:, 2), xyz_heads(:, 3)); hold on;
    scatter3(xyz_heads(:, 1), xyz_heads(:, 2), xyz_heads(:, 3), 50, 'filled'); hold on;    
    % legs
    for i = 1:2
        idxl = (i-1)*cnt_models+ii;
        xyz_legsf = Xf(idxl, :);
        xyz_legsf = reshape(xyz_legsf, [nbones_leg, 3]);
        xyz_legsf = [legsf{idxl}(1,:); xyz_legsf];
        plot3(xyz_legsf(:, 1), xyz_legsf(:, 2), xyz_legsf(:, 3)); hold on;
        scatter3(xyz_legsf(:, 1), xyz_legsf(:, 2), xyz_legsf(:, 3), 80, 'filled'); hold on;
    end
    for i = 1:2
        idxl = (i-1)*cnt_models+ii;
        xyz_legsf = Xb(idxl, :);
        xyz_legsf = reshape(xyz_legsf, [nbones_leg, 3]);
        xyz_legsf = [legsb{idxl}(1,:); xyz_legsf];
        plot3(xyz_legsf(:, 1), xyz_legsf(:, 2), xyz_legsf(:, 3)); hold on;
        scatter3(xyz_legsf(:, 1), xyz_legsf(:, 2), xyz_legsf(:, 3), 80, 'filled'); hold on;
    end
    % tail
    xyz_tail = reshape(Xt(ii,:), [nbones_t, 3]);
    xyz_tail = [tail{ii}(1,:); xyz_tail];
    plot3(xyz_tail(:, 1), xyz_tail(:, 2), xyz_tail(:, 3)); hold on;
    scatter3(xyz_tail(:, 1), xyz_tail(:, 2), xyz_tail(:, 3), 50, 'b', 'filled'); hold on;    

    % key points
    kpts = reshape(Xkpts(ii, :), [n_imptpts, 3]);
    plot3(kpts(:, 1), kpts(:, 2), kpts(:, 3)); hold on;
    scatter3(kpts(:, 1), kpts(:, 2), kpts(:, 3), 50, 'black', 'filled'); hold on;
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal;
    hold off;
end
%% write down
nkpts = n_imptpts + nbones_h + nbones_leg*4 + nbones_t;
xyz = zeros(nkpts,3);
for ii = 1:cnt_models
    % key points
    kpts = reshape(Xkpts(ii, :), [n_imptpts, 3]);
    % head
    xyz_heads = reshape(Xh(ii,:), [nbones_h, 3]);
    % legs
    xyz_legsf_l = reshape( Xf(ii, :), [nbones_leg, 3]);
    xyz_legsf_r = reshape( Xf(cnt_models+ii, :), [nbones_leg, 3]);
    xyz_legsb_l = reshape(Xb(ii, :), [nbones_leg, 3]);
    xyz_legsb_r = reshape(Xb(cnt_models+ii, :), [nbones_leg, 3]);
    % tail
    xyz_tail = reshape(Xt(ii,:), [nbones_t, 3]);

    xyz = [kpts; xyz_heads; xyz_legsf_l; xyz_legsf_r; xyz_legsb_l; xyz_legsb_r; xyz_tail];
    path_write = sprintf(formatSpec, 'horse_kpts',idx(ii));
    pcwrite(pointCloud(xyz), path_write);
end
