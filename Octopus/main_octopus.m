clear all;
close all;
%% �Ǽܺ���ȡ�Ǽܺ�õ��ľ���֮���
formatSpec = 'D:/GitHub/data/Octopus/skeletons/mesh_%.4d.ply';
W_formatSpec = 'D:/GitHub/data/Octopus/final/mesh_%.4d.ply';
P_formatSpec = 'D:/GitHub/data/Octopus/skeletons/mesh_%.4d_skeleton';
% Ŀǰֻ��5-40������, 56-58, 60:112
idx = [20];

% ��������
N = size(idx, 2);
cnt_leg = 8;

cnt_models = 0;
xyz = cell(N, 1);
xyz_models = cell(N,1);
parts_set = cell(N, cnt_leg);
parts_len = cell(N, cnt_leg);
for ii = idx
    paths = sprintf(formatSpec, ii);
    P_path = sprintf(P_formatSpec,  ii);
    ptCloud = pcread(paths);
      
    if ii< 47||ii==83
        xyz0 = ptCloud.Location;
        alphaY = -20./180*pi;
        transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
        tmp = transformY * xyz0';
        ptCloud = pointCloud(tmp');    
    end
    if ii == 86||ii==85
        xyz0 = ptCloud.Location;
        alphaY = 10./180*pi;
        transformY = [cos(alphaY), 0, sin(alphaY); 0, 1, 0; -sin(alphaY),  0, cos(alphaY)];
        tmp = transformY * xyz0';
        ptCloud = pointCloud(tmp');    
    end

    P = load(P_path);
    P = P.P;
    bbox = [ptCloud.XLimits; ptCloud.YLimits; ptCloud.ZLimits];
    rs = bbox(:,1)-bbox(:,2);
    diameter = sqrt(dot(rs,rs));

    % init key points
    [importantpts, parts, legslength]  = keypts_octopus(ptCloud.Location, P, 0);
    cnt_models = cnt_models+1;

%     xyz{cnt_models} = ptCloud.Location;
    xyz_model = ptCloud.Location;
    xyz_models{cnt_models} = xyz_model;
    % ֱ�Ӵ洢�ؼ��������
    imptpts{cnt_models} = xyz_model(importantpts, :);
    
    parts_set{cnt_models} = parts;
    parts_len{cnt_models} = legslength;
end

%ģ�����Ż���ʹ�ø�����һ������ ��̫�̵Ĳ���
nsample = 10;
nsample_head = 4;
parts_xyz = {cnt_models};
for m = 1:cnt_models
    %ÿ��ģ�͵ĸ���part�ĳ���
    lens_model = parts_len{m};
    parts_model = parts_set{m};
    max_lens = zeros(cnt_leg,1);
    xyz = xyz_models{m};%ģ�͵�����
    top_pts = zeros(cnt_leg,3);
    for i = 1:cnt_leg
        max_lens(i) = max(lens_model(i,:), [], 2);
        part_idx = parts_model{i};%ĳһ��part�ĵ���������Ӷ˵㿪ʼ
        ps = xyz(part_idx(1), :); %��һ������β��
        pe = xyz(part_idx(end), :);
        if dot(ps,ps) < dot(pe,pe)
            parts_model{i} = fliplr(part_idx);
        end
        part_idx = parts_model{i};
        top_pts(i,:) = xyz(part_idx(end), :);
    end
    sort_lens = sort(max_lens);
    avg_len = mean(sort_lens(end-6:end)); %ÿ���������ȡ�ĳ���
    middle_pts = mean(top_pts);
    parts_xyz_model = {cnt_leg};%������part��xyz�ļ���
    for i = 1:cnt_leg
        len_i = lens_model(i, :);%ĳһ��part�ĳ���
        len_i = len_i(len_i > 0);
        part_idx = parts_model{i};%ĳһ��part�ĵ���������Ӷ˵㿪ʼ
        part_idx = part_idx(1:end-3);
        len_i = len_i(1:end-1);
        if (0) %(max(len_i) < avg_len) %�������Ĳ�������һ��������֤����
%             len_tmp = max(size(part_idx,2)-4,3);
            pts = xyz(part_idx, :);
            [coefs,~,latent] = pca(pts);
            if (size(coefs,2) < 1)
                part_idx
            end
            dir = coefs(:, 1)';%������ 3*1
            dir = dir/norm(dir);
            if dot((pts(end-1,:)-pts(end,:)), dir) > 0
                dir = dir*-1;
            end
            pt_last = xyz(part_idx(end), :);
            newpt = dir*(avg_len-max(len_i)) + pt_last;
            parts_xyz_model{i} = [xyz(part_idx, :); newpt];
            len_i = [len_i,avg_len];
        else
            parts_xyz_model{i} = xyz(part_idx, :); 
        end
        
        parts_xyz_model{i} = [parts_xyz_model{i};middle_pts];
        %%����������ɵȳ��ĳ���
        parts_xyz_model{i} = resampleByLength(parts_xyz_model{i}, nsample);
    end
    head_idx = parts_model{9};
    idx_end = max(size(head_idx,2)-2,1);
    head_idx = head_idx(1:idx_end);
    parts_xyz_model{9} = [xyz(head_idx,:); middle_pts];
    parts_xyz_model{9} = resampleByLength(parts_xyz_model{9}, nsample_head);
    parts_xyz{m} = parts_xyz_model;
    
end

% % ��ʾ�ز������ģ��
% for ii = 1:cnt_models
%     figure(idx(ii));
%     parts = parts_xyz{ii};
%     for j = 1:cnt_parts
%         p = parts{j};
%         %scatter3(p(:,1), p(:,2), p(:,3), 15, 'filled'); hold on;
%         plot3(p(:,1), p(:,2), p(:,3), '.-'); hold on;
%     end
%     
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     axis equal;
%     hold off;
% end

%% ��������ȶ��룬�����е�ģ�Ͷ�������ͬ��˳������
npart = size(parts_set, 1);
aligned_xyz_models = {};
for i = 1:cnt_models
    cell_xyz = parts_xyz{i};
	legs = align_legs_octopus( cell_xyz(1:8), cnt_leg);
    legs{9} = cell_xyz{9}; %���һ����ͷ
	aligned_xyz_models{i} = legs;
end

%% plot
C = jet(cnt_leg+1);
%C = hsv(cnt_leg+1);
for ii = 1:cnt_models
    h1 = figure(idx(ii)+100);
    set( h1, 'Color','w')
    aligned_xyz = aligned_xyz_models{ii}; 
%     scatter3(xyz_model(:,1), xyz_model(:,2), xyz_model(:,3), 15, 'filled'); hold on;

    for j = 1:cnt_leg+1
        p = aligned_xyz{j};
        %scatter3(p(:,1), p(:,3), p(:,2), 15, 'filled'); hold on;
        plot3(p(:,3), p(:,1), p(:,2), '.-', 'LineWidth', 5, 'MarkerSize', 42, 'Color', C(j,:)); hold on;
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal;
    axis off;
    hold off;
end
% %% write down
% n_pts = nsample*cnt_leg+nsample_head;
% for ii = 1:cnt_models
%     % key points
%     pts = zeros(n_pts, 3);
%     aligned_xyz = aligned_xyz_models{ii}; 
%     for j = 1:cnt_leg
%         p = aligned_xyz{j};
%         n1 = j*nsample-nsample+1;
%         n2 = j*nsample;
%         pts(n1:n2,:) = p;
%     end
%     pts(end-nsample_head+1:end,:) = aligned_xyz{9};
%     path_write = sprintf(W_formatSpec, idx(ii));
%     pcwrite(pointCloud(pts), path_write);
% end

% %% write edge down
% n_pts = nsample*cnt_leg+nsample_head;
% 
% N = zeros(n_pts, n_pts);
%     for j = 1:cnt_leg
%         n1 = j*nsample-nsample+1;
%         n2 = j*nsample;
%         for i = n1:n2-1
%             N(i,i+1) = 1;
%         end
%     end
% for i = n_pts-nsample_head+1:n_pts-1
% 	N(i,i+1) = 1;
% end
% [X,Y] = find(N);
% edges = [X';Y']'
% path_edge_write = 'D:/GitHub/data/Octopus/final/edges.txt';
% fp=fopen(path_edge_write,'a');%'A.txt'Ϊ�ļ�����'a'Ϊ�򿪷�ʽ���ڴ򿪵��ļ�ĩ��������ݣ����ļ��������򴴽���
% for e  = 1:size(edges,1)
%     fprintf(fp,'%d %d\n',edges(e,1), edges(e,2));%fpΪ�ļ������ָ��Ҫд�����ݵ��ļ���ע�⣺%d���пո�
% end
% 
% fclose(fp);%�ر��ļ���


