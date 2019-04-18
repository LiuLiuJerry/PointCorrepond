function [importantpts, parts, partslength] = keypts_horse(xyz, P, bplot)
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
% % 最小生成树，找边界点
    dists = pdist2(xyz, xyz);
%     gra = graph(dists); %距离作为边的权重
%     G = minspantree(gra);% 最小生成树
    self =  ~diag(ones(N,1));
    
    A = P.spls_adj & self;
    Gra = graph(A);
    D = degree(Gra);
    
    % 删除多余的边
    neighbors = {};
    for i = 1:N       
        idx = find(A(:,i));
        dis = dists(idx, i);
        nbs = [double(idx), dis]; %idx是列向量
        nbs = sortrows(nbs,2); % 第一遍排序是为了保留最近的点
        
        %test: 将从中心指向邻居的向量给降维
        if  D(i) < 2
            neighbors{i} = nbs;
            continue;
        end
        diri = xyz(i,:)-xyz(nbs(:,1),:);
        for j = 1:size(nbs,1) %计算每个走向的方向
            diri(j, :) = diri(j, :)./norm(diri(j,:));              
        end
        for j = 2:size(nbs, 1) %第二个点的走向应该跟前面几个点不一样
            for k = 1:j-1
                if dot(diri(j,:), diri(k,:)) > 0.98
                    idxj = (nbs(j,1));
                    idxk = (nbs(k,1));
                    A(idxj,i) = 0;
                    A(i,idxj) = 0;
                end
            end
        end
        idx = find(A(:,i));
        dis = dists(idx, i);
        nbs = [double(idx), dis]; %idx是列向量
        nbs = sortrows(nbs,2);
        
        neighbors{i} = nbs;
        
%         if size(nbs,1) > 2
%             importantpts = [importantpts; i];
%         end
    end
    D = sum(A, 2);   
    
    %% 初步按拓扑分割
    idx_Start = find(D<2);
    npart = size(idx_Start, 1);
    parts = cell(npart+1,1);
    partslength = zeros(npart,30);
    seg = zeros(size(xyz,1),1);
    for t = 1:npart
        idx = idx_Start(t);
        pre = 0;
        i = 1;
        part = zeros(1, 20);
        while D(idx) < 3
            part(i) = idx;
            seg(idx) = t;
            nei = neighbors{idx};
            if D(idx) == 1 && nei(1,1) == pre
                break;
            end
            if nei(1,1)~=pre
                pre = idx;
                idx = nei(1,1);
                d = max(1, i-1);
                partslength(t, i) = partslength(t, d) + nei(1,2);
            else 
                pre = idx;
                idx = nei(2,1);
                d = max(1, i-1);
                partslength(t, i) = partslength(t, d) + nei(2,2);
            end
            i = i+1;
        end
        part(i) = idx; %每个部分的最右一个点都是连接点
        seg(idx) = t;
        part = part(part > 0);
        parts{t} = part;
    end
    
    %%
    %删掉一些部分，合并一些部分
    pos = [[1:npart]', xyz(idx_Start,:), max(partslength,[],2)];
    pos = sortrows(pos, 5); 
    maxlen = max(pos(:,5));
    details = pos(pos(:,5) < maxlen * 0.4, :);
    details = sortrows(details, 3);%y最大的是头
    idx_head = details(end,1);
    details = details(1:end-1,:);
    nde = size(details, 1);
    delta_y = max(xyz(idx_Start,2))-min(xyz(idx_Start,2));
    py = min(xyz(idx_Start,2)) + delta_y*0.2;
    for i = 1:nde %删除其中的一些部分
        idxp = details(i,1);
        idxpts = parts{idxp};
        pts_d = xyz(idxpts,:);
        %if s_y < 0: %是脚上的问题
        dir = pts_d(1,:)-pts_d(end,:);
        dir = dir./norm(dir);
        cos_d1 = dot(dir,[0,-1,0]);   
        cos_d2 = dot(dir,[0,0,-1]);
        if ( cos_d2>0) || (pts_d(end,2)>py)%是没有意义的点，需要删除
            idxnot = idxpts(end);
            idxnot_nei = idxpts(end-1);
            A(idxnot,idxnot_nei) = 0;
            A(idxnot_nei,idxnot) = 0;
            nei = neighbors{idxnot};
            nei = nei(nei(:,1) ~= idxnot_nei,:);
            neighbors{idxnot}= nei;
            nei = neighbors{idxnot_nei};
            nei = nei(nei(:,1) ~= idxnot,:);
            neighbors{idxnot_nei}= nei;            
            
            details(i,:)=0;
            %idx_Start(idxp) = 0;
            continue;
        end
    end
    details = details(details(:,1)>0,:);
    D = sum(A, 2);
    
%     nde = size(details, 1);
%     %idx_Start = idx_Start(idx_Start>0);
%     for j = 1:nde %剩下的部分重新开始连接
%         idxp = details(j,1);
%         idxpts = parts{idxp};
%         
%         idxnot = idxpts(end);
%         idxnot_nei = idxpts(end-1);
%         %重新开始计算
%         idx = idxnot;
%         pre = idxnot_nei;
%         i = size(idxpts, 2)+1;
%         part = parts{idxp};
%         t = idxp;
%         while D(idx) < 3
%             part(i) = idx;
%             seg(idx) = t;
%             nei = neighbors{idx};
%             if D(idx) == 1 && nei(1,1) == pre
%                 break;
%             end
%             if nei(1,1)~=pre
%                 pre = idx;
%                 idx = nei(1,1);
%                 d = max(1, i-1);
%                 partslength(t, i) = partslength(t, d) + nei(1,2);
%             else 
%                 pre = idx;
%                 idx = nei(2,1);
%                 d = max(1, i-1);
%                 partslength(t, i) = partslength(t, d) + nei(2,2);
%             end
%             i = i+1;
%         end
%         part(i) = idx; %每个部分的最右一个点都是连接点
%         seg(idx) = t;
%         part = part(part > 0);
%         parts{idxp} = part;
%     end

    %% 分割成头尾和四条腿
    idx_Start = find(D<2);
    npart = size(idx_Start, 1);
    parts = cell(npart+1,1);
    partslength = zeros(npart,30);
    seg = zeros(size(xyz,1),1);
    for t = 1:npart
        idx = idx_Start(t);
        pre = 0;
        i = 1;
        part = zeros(1, 20);
        while D(idx) < 3
            part(i) = idx;
            seg(idx) = t;
            nei = neighbors{idx};
            if size(nei,1) == 0|| (D(idx) == 1 && nei(1,1) == pre)
                break;
            end
            if nei(1,1)~=pre
                pre = idx;
                idx = nei(1,1);
                d = max(1, i-1);
                partslength(t, i) = partslength(t, d) + nei(1,2);
            else 
                pre = idx;
                idx = nei(2,1);
                d = max(1, i-1);
                partslength(t, i) = partslength(t, d) + nei(2,2);
            end
            i = i+1;
        end
        part(i) = idx; %每个部分的最右一个点都是连接点
        seg(idx) = t;
        part = part(part > 0);
        parts{t} = part;
    end
    %%    
	pos = [[1:npart]', xyz(idx_Start,:), max(partslength,[],2)];
    pos = sortrows(pos, 5); 
    
    maxlen = max(pos(:,5));
    pos = sortrows(pos, 5); 
%     phead = pos(pos(:,5) < maxlen * 0.4, :);
%     if isempty(phead)
%         idx_head = [];
%     else
%         idx_head = phead(end,1);
%     end

    fpos = pos(end-3:end,:);
    pos = sortrows(pos, 4);
    
     %先 z排序 再x排序 x-z:(-1,-1),(-1,1),(1,-1),(1,1)右后，右前，左后，左前   
    fpos = sortrows(fpos, 3); % x排序 结果为从第一象限开始逆时针计算
    fpos(1:2, :) = sortrows(fpos(1:2,:), 2); %左胳膊右胳膊
    fpos(3:4, :) = sortrows(fpos(3:4,:), 2);%左腿右腿

    idx_legs = fpos(:,1);
    
%     idx_p = [idx_ear;idx_head; idx_legs(:); idx_tail];
    idx_p = [idx_head; idx_legs(:)];
    parts = parts(idx_p);
    % 剩余的是脊椎骨架
    backbone_idx = find(seg == 0) ;
    node_arml = parts{2};
    node_armr = parts{3};
    node_legl = parts{4};
    node_legr = parts{5};
    parts{end+1} = [node_arml(end),node_armr(end), reshape(backbone_idx, [1, size(backbone_idx,1)]), node_legl(end), node_legr(end)];
    parts{end} = unique(parts{end});
    %parts{end+1} = reshape(backbone_idx, [1, size(backbone_idx,1)]);
    % 关键点
    importantpts = find(D > 2);
    ipos = [importantpts, xyz(importantpts, :)];
    ipos = sortrows(ipos, 3); 
    ipos = sortrows(ipos, 4); %从小到大
    importantpts = ipos(:,1);

    %%
npart = size(parts,1);
if bplot
    figure()
%     scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 15, 'filled'); hold on;
    plot_skeleton(xyz, A);hold on;
    for i = 1:npart
        idx = parts{i};
        pts = xyz(idx, :);
        scatter3(pts(:,1), pts(:,2), pts(:,3),  48, 'filled'); hold on;
    end
     idx_Start = xyz(D<2,:);
%     D_3 = xyz(D>2,:);
%      scatter3(idx_Start(:,1), idx_Start(:,2), idx_Start(:,3), 72, 'filled'); hold on;
%     scatter3(D_3(:,1), D_3(:,2), D_3(:,3), 72, 'b', 'filled'); hold on;
    axis equal;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off;
end

end
            
                
    
    
    
    

