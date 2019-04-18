clear all;
close all;
%% 骨架和提取骨架后得到的矩阵啊之类的
formatSpec = 'D:/GitHub/data/human/I_squat/skeletons/mesh_%.4d.ply';
P_formatSpec = 'D:/GitHub/data/human/I_squat/skeletons/mesh_%.4d_skeleton';
out_formatSpec = 'D:/GitHub/data/human/I_squat/squat_kpts/mesh_%.4d.ply';
idx = [0:49];

% 定义骨架的结构
n_imptpts = 3;
nbones_h = 1;
nbones_leg = 3;
nbones_arm = 2;

% 声明数据
N = size(idx, 2);
imptpts = cell(N, 1);
ear = cell(N,1);
head = cell(N, 1);
arms_r = cell(N, 1);
legs_r = cell(N, 1);
arms_l = cell(N, 1);
legs_l = cell(N, 1);
backbone = cell(N,1);

cnt_models = 0;
xyz = cell(N, 1);
xyz_models = cell(N,1);
for ii = idx
    paths = sprintf(formatSpec, ii);
    P_path = sprintf(P_formatSpec, ii);
    ptCloud = pcread(paths);
    P = load(P_path);
    P = P.P;
    bbox = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
    rs = bbox(:,1)-bbox(:,2);
    diameter = sqrt(dot(rs,rs));

    % init key points
    [importantpts, parts, legslength]  = keypts_human(ptCloud.Location, P, 1);
    cnt_models = cnt_models+1;

    xyz{cnt_models} = ptCloud.Location;
    xyz_model = ptCloud.Location;
    xyz_models{cnt_models} = xyz_model;
    % 直接存储坐标
    imptpts{cnt_models} = xyz_model(importantpts, :);%z轴从小到大排序
    % ear{cnt_models} = xyz_model(fliplr(parts{1}), :);
    head{cnt_models} = xyz_model(fliplr(parts{1}), :);
    % (-1,-1),(1,-1),(-1,1),(1,1)右后，左后，右前，左前
    arms_l{cnt_models} = xyz_model(fliplr(parts{4}), :);
    arms_r{cnt_models} = xyz_model(fliplr(parts{5}), :);
    legs_l{cnt_models} = xyz_model(fliplr(parts{2}), :);
    legs_r{cnt_models} = xyz_model(fliplr(parts{3}), :); %反转向量，从里往外
    
    backbone{cnt_models} = xyz_model(fliplr(parts{6}), :);
end

%% 优化关键点
% 用开始提取的关键点做x0
x0_kpts = zeros(cnt_models, n_imptpts*3);
for i = 1:cnt_models
    x0_tmp = imptpts{i}([1,floor(end/2),end]', :);
    x0_kpts(i, :) = x0_tmp(:);
end
% 先优化脊椎骨的关节点
Xkpts = fmincon(@(X)fun_impdist(X, backbone, imptpts, n_imptpts), x0_kpts, [], [], [], [], [], [], @(X)nonlcon_imp(X,  imptpts, n_imptpts));

%% 优化骨节连接点
% 估算x0
xh0 = zeros(cnt_models, nbones_h*3);
xf0 = zeros(cnt_models*2, nbones_arm*3);
xb0 = zeros(cnt_models*2, nbones_leg*3);
for t = 1:cnt_models %把腿和脊椎的节点组织起来
    kpts = reshape(Xkpts(t, :), [n_imptpts, 3]);
    kpts = sortrows(kpts,3);
    head{t} = [kpts(end,:); head{t}];  %这时候是按照从里到外的顺序组织的，包括后面存储的点也是
    arms_l{t} = [kpts(end,:);arms_l{t}];
    arms_r{t} = [kpts(end,:);arms_r{t}];
    legs_l{t} = [kpts(1,:);legs_l{t}];
    legs_r{t} = [kpts(1,:);legs_r{t}];    
end
    
for t = 1:cnt_models %初始化这些值
    init_h = head{t}(end, :);
    xh0(t, :) = init_h(:);
    initarm_r = arms_r{t}([floor(end/2), end]', :);
    initarm_l = arms_l{t}([floor(end/2), end]', :);
    xf0(t,:) = initarm_l(:); %先排n个左腿，再排n个右腿
    xf0(cnt_models+t, :) = initarm_r(:);
    initleg_r = legs_r{t}([floor(end/2), end-1:end]', :);
    initleg_l = legs_l{t}([floor(end/2), end-1:end]', :);
    xb0(t,:) = initleg_l(:);
    xb0(cnt_models+t, :) = initleg_r(:);
end

VLB = [];
VUB = [];
% 应该把关节点也加入坐标中
arms = [arms_l, arms_r];
legsb = [legs_l, legs_r];
% 传入坐标以前腿的坐标应该减去最开始的坐标
% 优化腿和头的关节点
Xh = fmincon(@(X)fun_dist(X, head, nbones_h), xh0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, head, nbones_h));
Xf = fmincon(@(X)fun_dist(X, arms, nbones_arm), xf0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, arms, nbones_arm));
Xb = fmincon(@(X)fun_dist(X, legsb, nbones_leg), xb0, [], [], [], [], VLB, VUB, @(X)nonlcon(X, legsb, nbones_leg));
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
        xyz_arms = Xf(idxl, :);
        xyz_arms = reshape(xyz_arms, [nbones_arm, 3]);
        xyz_arms = [arms{idxl}(1,:); xyz_arms];
        plot3(xyz_arms(:, 1), xyz_arms(:, 2), xyz_arms(:, 3)); hold on;
        scatter3(xyz_arms(:, 1), xyz_arms(:, 2), xyz_arms(:, 3), 80, 'filled'); hold on;
    end
    for i = 1:2
        idxl = (i-1)*cnt_models+ii;
        xyz_legsf = Xb(idxl, :);
        xyz_legsf = reshape(xyz_legsf, [nbones_leg, 3]);
        xyz_legsf = [legsb{idxl}(1,:); xyz_legsf];
        plot3(xyz_legsf(:, 1), xyz_legsf(:, 2), xyz_legsf(:, 3)); hold on;
        scatter3(xyz_legsf(:, 1), xyz_legsf(:, 2), xyz_legsf(:, 3), 80, 'filled'); hold on;
    end

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
nkpts = n_imptpts-1 + nbones_h + nbones_leg*2 + nbones_arm*2 ;
xyz = zeros(nkpts,3);
for ii = 1:cnt_models
    % key points
    kpts = reshape(Xkpts(ii, :), [n_imptpts, 3]);
    kpts = kpts(end-1:end,:);
    % head
    xyz_heads = reshape(Xh(ii,:), [nbones_h, 3]);
    % legs
    xyz_legsf_l = reshape( Xf(ii, :), [nbones_arm, 3]);
    xyz_legsf_r = reshape( Xf(cnt_models+ii, :), [nbones_arm, 3]);
    xyz_legsb_l = reshape(Xb(ii, :), [nbones_leg, 3]);
    xyz_legsb_r = reshape(Xb(cnt_models+ii, :), [nbones_leg, 3]);

    xyz = [kpts; xyz_heads; xyz_legsf_l; xyz_legsf_r; xyz_legsb_l; xyz_legsb_r; ];
    path_write = sprintf(out_formatSpec, idx(ii));
    pcwrite(pointCloud(xyz), path_write);
end
