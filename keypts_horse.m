function [importantpts, parts, partslength] = keypts_horse(xyz, P, bplot)
%% parameters

% thre: thresold for important points
% r: search redius

% %���ļ�
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
% % ��С���������ұ߽��
    dists = pdist2(xyz, xyz);
%     gra = graph(dists); %������Ϊ�ߵ�Ȩ��
%     G = minspantree(gra);% ��С������
    self =  ~diag(ones(N,1));
    
    A = P.spls_adj & self;
    Gra = graph(A);
    D = degree(Gra);
    
    % ɾ������ı�
    neighbors = {};
    for i = 1:N       
        idx = find(A(:,i));
        dis = dists(idx, i);
        nbs = [double(idx), dis]; %idx��������
        nbs = sortrows(nbs,2); % ��һ��������Ϊ�˱�������ĵ�
        
        %test: ��������ָ���ھӵ���������ά
        if  D(i) < 2
            neighbors{i} = nbs;
            continue;
        end
        diri = xyz(i,:)-xyz(nbs(:,1),:);
        for j = 1:size(nbs,1)
            diri(j, :) = diri(j, :)./norm(diri(j,:));              
        end
        for j = 2:size(nbs, 1)
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
        nbs = [double(idx), dis]; %idx��������
        nbs = sortrows(nbs,2);
        
        neighbors{i} = nbs;
        
%         if size(nbs,1) > 2
%             importantpts = [importantpts; i];
%         end
    end
    D = sum(A, 2);   
    
    %% �ָ��ͷβ��������
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
        part(i) = idx;
        seg(idx) = t;
        part = part(part > 0);
        parts{t} = part;
    end
    pos = [[1:npart]', xyz(idx_Start,:), max(partslength,[],2)];
    maxlen = max(pos(:,5));
    pos = sortrows(pos, 5);
    phead = pos(pos(:,5) < maxlen * 0.4, :);
    if isempty(phead)
        idx_head = [];
    else
        idx_head = phead(end,1);
    end
    pos = pos(pos(:,5) > maxlen * 0.4, :);
    
    pos = sortrows(pos, 4);
    idx_tail = pos(1,1);
    fpos = pos(2:5,:);
     %�� z���� ��x���� x-z:(-1,-1),(-1,1),(1,-1),(1,1)�Һ���ǰ�������ǰ   
    fpos = sortrows(fpos, 2); % x���� ���Ϊ�ӵ�һ���޿�ʼ��ʱ�����
    fpos(1:2, :) = sortrows(fpos(1:2,:), 4);
    fpos(3:4, :) = sortrows(fpos(3:4,:), 4);

    idx_legs = fpos(:,1);
    
%     idx_p = [idx_ear;idx_head; idx_legs(:); idx_tail];
    idx_p = [idx_head; idx_legs(:); idx_tail];
    parts = parts(idx_p);
    % ʣ����Ǽ�׵�Ǽ�
    parts{npart+1} = find(seg == 0);
    % �ؼ���
    importantpts = find(D > 2);
    ipos = [importantpts, xyz(importantpts, :)];
    ipos = sortrows(ipos, 3); 
    ipos = sortrows(ipos, 4); %��С����
    importantpts = ipos(:,1);

    %%
if bplot
    figure()
%     scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 15, 'filled'); hold on;
    plot_skeleton(xyz, A);
    for i = 1:npart+1
        idx = parts{i};
        pts = xyz(idx, :);
        scatter3(pts(:,1), pts(:,2), pts(:,3),  48, 'filled'); hold on;
    end
    idx_Start = xyz(D<2,:);
    D_3 = xyz(D>2,:);
    scatter3(idx_Start(:,1), idx_Start(:,2), idx_Start(:,3), 72, 'filled'); hold on;
    scatter3(D_3(:,1), D_3(:,2), D_3(:,3), 72, 'b', 'filled'); hold on;
    axis equal;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off;
end

end
            
                
    
    
    
    

