function [importantpts, parts, partslength] = keypts_octopus(xyz, P, bplot)
%% parameters
    %importantpts �� ���еĹǼ��ϵĵ���������һ��
    %parts �� �ָ�ɸ������ֵ�����


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
    
    % ���ھӣ� ɾ������ı�
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
    
    %% �ָ�n����
    idx_Start = find(D<2); %��ʼ��ֻ��һ���ھ�
    npart = size(idx_Start, 1);
    parts = cell(npart+1,1);
    partslength = zeros(npart,30);
    seg = zeros(size(xyz,1),1); %ÿ���������ĸ�part�ı��
    for t = 1:npart
        idx = idx_Start(t);
        pre = 0;
        i = 1;
        part = zeros(1, 20); %��¼ÿ��part ������
        while D(idx) < 3%��������
%             if i > 2 %���ǰ���ĽǶ����̫���break��
%                 v1 = part(i-1)-part(i-2);
%                 v2 = part(i)-part(i-1);
%                 ang = dot(v1,v2)/(norm(v1)*norm(v2));
%                 if(ang < cos(60/180*pi) )
%                     break;
%                 end
%             end
            part(i) = idx;
            seg(idx) = t;
            nei = neighbors{idx};
            if D(idx) == 1 && nei(1,1) == pre %�ߵ�ͷ��
                break;
            end
            if nei(1,1)~=pre %û���߹����ھ�
                pre = idx;
                idx = nei(1,1);
                d = max(1, i-1);
                partslength(t, i) = partslength(t, d) + nei(1,2);%���³���
            else 
                pre = idx;
                idx = nei(2,1);
                d = max(1, i-1);
                partslength(t, i) = partslength(t, d) + nei(2,2);
            end
            i = i+1;
        end 
        
        if seg(idx) == 0
            part(i) = idx;%���һ���ҵ����ھ��Ƿֲ��
            seg(idx) = t;
        else
            part = part(part > 0);
            if size(part,2) < 5
                partslength(t, i-1) = 0;
            else
                tmp = part(1,end-3:end); %ȥ�������Сβ��
                pts = xyz(tmp', :);
                [coefs,~,latent] = pca(pts);
                dir = coefs(:, 1)';%������ 3*1
                main_dir = dir/norm(dir);
                if pts(end,:)*main_dir'<0 %���һ������ͶӰӦ��������
                   main_dir = -main_dir;
                end 
                 v1 = xyz(idx,:) - xyz(part(1,end),:);
                 v1 = v1/norm(v1);
                if dot(main_dir,v1) > cos(1/3*pi) %ƫ�����������60�㣬˵�����԰ѷֲ����ɽ���
                     part(i) = idx;%���һ���ҵ����ھ��Ƿֲ��
                     seg(idx) = t;
                else
                     partslength(t, i-1) = 0;
                end
            end
                      
        end
        
%             if size(part,2)>6 %���һ���ֵ�������
%                 tmp = part(1,end-5:end); %ȥ�������Сβ��
%                 pts = xyz(tmp', :);
%                 [coefs,~,latent] = pca(pts);
%                 dir = coefs(:, 1)';%������ 3*1
%                 main_dir = dir/norm(dir);
%                 if pts*main_dir'<0
%                     main_dir = -main_dir;
%                 end 
%             end

        part = part(part > 0);
        parts{t} = part;
    end
    % ÿ��part�ı�ţ� x, y, z, part����, ���ڶ�ά�����ֵ
    pos = [[1:npart]', xyz(idx_Start,:), max(partslength,[],2)];
    pos = flipud(sortrows(pos, 5));%���ݳ�������
 
    npart = 9;
    idx_p = pos(1:npart,1);
    all_parts =  parts;
    %���
    parts = parts(idx_p);
    partslength = partslength(idx_p, :);
    
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
%     npart = size(parts, 1);
%     plot_skeleton(xyz, A);
    for i = 1:npart
        idx = parts{i};
        pts = xyz(idx, :);
        scatter3(pts(:,1), pts(:,2), pts(:,3),  48, 'filled'); hold on;
    end
%     idx_Start = xyz(D<2,:);
%     D_3 = xyz(D>2,:);
%     scatter3(idx_Start(:,1), idx_Start(:,2), idx_Start(:,3), 72, 'filled'); hold on;
%     scatter3(D_3(:,1), D_3(:,2), D_3(:,3), 72, 'b', 'filled'); hold on;
    axis equal;
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold off;
end

end
            
                
    
    
    
    

