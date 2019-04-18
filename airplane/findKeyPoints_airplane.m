function [xyz, seg, dir, parts, importantpts, keypts_idx, boundary_idx] = findKeyPoints(xyz, P, thre, r)
%% parameters

% thre: thresold for important points
% r: search redius

% %读文件
% clear all;
% close all;
% formatSpec = 'D:\\GitHub\\data\\SkeletonsbyJerry\\%s\\airplane_0%d.ply';
% filename = 'airplane';
% filename_ort = 'airplane_ort';
% filename_ali = 'airplane_aligned';

% path = sprintf(formatSpec, filename_ort, ii);
% name = sprintf('airplane_0%d,ply', ii);
% ptCloud = pcread(path);
% xyz = ptCloud.Location;

    N = size(xyz,1);
% 最小生成树，找边界点
    dists = pdist2(xyz, xyz);
    gra = graph(dists); %距离作为边的权重
    G = minspantree(gra);% 最小生成树
    A = P.spls_adj;
    D = degree(G);
    
    %计算点的方向 pca
    importantpts = [];
    % 可以换为serchradius
    k = 5;
    [idx_n, ~] = knnsearch(xyz, xyz, 'K', k);
    nei_idx = (dists < r);
%     nei_idx = nei_idx | A;
    nei_idx = A | diag(ones(N, 1));

    dir = zeros(N, 3);
    neighbors = {};
%     % 找邻居
%     for i = 1:N
%         idx = find(nei_idx(:,i));
%          %度为1才可以至少有一个邻居，否则至少有两个邻居
%          if size(idx, 1) < 3 && D(i) > 1
%              idx = idx_n(i, 1:3);
%              nei_idx(i, idx(2)) = 1;
%              nei_idx(idx(2), i) = 1;   
%              nei_idx(i, idx(3)) = 1;
%              nei_idx(idx(3), i) = 1; 
%          elseif size(idx,1) < 2 && D(i) == 1
%             idx = idx_n(i, 1:2);
%             nei_idx(i, idx(2)) = 1;
%             nei_idx(idx(2), i) = 1;
%   
%         end       
%     end
    for i = 1:N       
        % 最少取一个邻居就可以了，太远的不要
%         [idx, dis] = findNeighborsInRadius(pointCloud(xyz), xyz(i,:), r, 'Sort', true);
        idx = find(nei_idx(:,i));
        dis = dists(idx, i);
        nbs = [double(idx), dis]; %idx是列向量
        nbs = sortrows(nbs,2);

        if size(idx,1) <= 2
            D(i) = 1;
        end

        neighbors{i} = nbs(2:end, :);
        neighbor = xyz(idx(:), :);
        %test: 将从中心指向邻居的向量给降维
        [coefs,~,latent] = pca(neighbor);
%         v = xyz(i)-neighbor;
%         dir(i, :) = mean(v,1);
        dir(i, :) = coefs(:, 1)';
        sum_la = sum(latent);
        if (latent(1)/sum_la) < thre
            importantpts = [importantpts; i];
%             key_nei = [key_nei; idx];
        end

    end

    %邻居数目为1的关键点消除
    importantpts = setdiff(importantpts, find(D<2));    

%%展示方向的计算结果
    % plot
    figure();
    for jj=1:N
        pts = [xyz(jj,:); xyz(jj,:)+dir(jj, :)./40];
        plot3(xyz(jj,1), xyz(jj,2), xyz(jj,3), 'o', 'MarkerSize', 3);
        plot3(pts(:,1), pts(:,2), pts(:,3)); hold on;
    end
    axis equal;
    title('Direction');
    hold off;

%法向或者曲率骤变的地方进行分割

%     idx_S = key_nei(:); %从关键点的邻居开始
%     idx_S = setdiff(idx_S, keypts);
    t = 1;
    seg = zeros(N,1);
    seg(importantpts,:) = -1;
    for i = 1:N
        dir_n(i, :) = dir(i, :)./norm(dir(i,:));
    end
%深搜
    %关键点聚类取平均,或者邻接的关键点合并成一个
    key_xyz = xyz(importantpts, :);
    count = zeros(size(importantpts));%列向量
    iptpts_count = zeros(size(importantpts));%列向量
    for i = 1:size(importantpts,1)
        tmp = neighbors{importantpts(i)};
        tmp = tmp(:,1);
        count(i) = size(tmp,1);
        iptpts_count(i) = size( setdiff(tmp, tmp(seg(tmp)>0)), 1);
    end
%     T = clusterdata(key_xyz, 1.2);
%     if size(unique(T),1) < 2
    init_K = iptpts_count(iptpts_count > 5);
    T = kmeans(key_xyz, 4, 'Start', 'uniform','Replicates', 5);
%     end
    parts = [];
    for ii = importantpts'
        key_nei = neighbors{ii};
        for jj = key_nei(:,1)' %关键点的邻居，不一定是关键点
            if seg(jj) > 0 
                continue; 
            end
            if seg(jj) == -1 && (T(importantpts == jj) == T(importantpts==ii))
                continue; % 标记在一起的关键点会被聚类到一起
            end
            seg(ii) = t;
            init_nei = neighbors{jj};
            idx_nei = init_nei(:,1);
            tmp = find(seg(idx_nei) == 0);
            nbi = init_nei(tmp, :);
            nbs = neighbors;
            nbs{jj} = nbi;
            [seg, pend] = tag(seg, dir_n, nbs, jj, t);       
            t = t+1;
            %保存一次搜索的两个顶点
            if pend ~= ii
                parts = [parts; [ii, pend]];
            end
        end
    end
    n_tags = t-1;
    notfound = xyz(seg == 0, :);

    size_T = size(unique(T), 1);
    key_points = zeros(size_T, 3);;
    keypts_idx = zeros(size_T,1);
    for t_k = 1:size_T
        idx_t = find(T == t_k); %关键点中的序号
        key_points(t_k, :) = mean(key_xyz(idx_t, :), 1);
        % 更新直线的两个端点
        keyidx_T = importantpts(idx_t);
        tmp = ismember(parts, keyidx_T);
        parts(tmp) = N + t_k; 
        keypts_idx(t_k) = N + t_k;
    end
    parts = unique(parts, 'rows');
    xyz = [xyz; key_points];
    dir = [dir; zeros(size(key_points,1),3)];
    boundary_idx = find(D<2);
    
    %之前的关键点合并到临近的segment里面去
    untagedkeys = [importantpts, count];
    untagedkeys = sortrows(untagedkeys, 2);
    untagedkeys = untagedkeys(:,1);
    for i = untagedkeys'
        nei = neighbors{i}; nei = nei(:,1);
        t = seg(nei);
        nei = nei(t > 0);
        un_nei = setdiff(nei, untagedkeys);
        if ~isempty(un_nei)
            seg(i) = seg(un_nei(1)); 
        elseif  ~isempty(nei)
            seg(i) = seg(nei(1)); 
        end
    end  
    seg = [seg; ones(size_T,1).*-1];
    
    % plot
    figure();
    scatter3(key_xyz(:, 1), key_xyz(:, 2), key_xyz(:, 3),20, T,'filled');  hold on;
    scatter3(key_points(:, 1), key_points(:, 2), key_points(:, 3),80,'filled');  hold on;
    if ~isempty(notfound)
    scatter3(notfound(:, 1), notfound(:, 2), notfound(:, 3),80,'filled', 'r');  hold on;
    end
%     scatter3(xyz(boundary_idx, 1), xyz(boundary_idx, 2), xyz(boundary_idx, 3), 80,'filled', 'o'); hold on;
    for ii = 1:n_tags
        seg_xyz = xyz(seg==ii, :);
        plot3(seg_xyz(:, 1), seg_xyz(:, 2), seg_xyz(:, 3), '.', 'MarkerSize', 6); hold on;
    end
    axis equal;
%     title('Segment and key points');
%     hold off;
%     
%     figure();
    for ii = parts'
        plot3(xyz(ii,1), xyz(ii,2), xyz(ii,3)); hold on;
    end
    axis equal;
    title('extracted lines');
    hold off;
end