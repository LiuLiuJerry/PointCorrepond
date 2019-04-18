function aligned_xyz = align_legs_octopus(xyz, cnt_parts)

aligned_idx = zeros(1,cnt_parts);
 %�ҵ���һ���ȣ����Ϊp
max_dot = 0;
p = 1;
init = [0,0,1];
 % ��z�� ������
for i = 1:cnt_parts
    pmean = mean(xyz{i});
    pmean = pmean./norm(pmean);
    if dot(pmean,init) > max_dot & pmean(1)<0
        max_dot = pmean(:,3);
        p = i;
    end
end
aligned_idx(1) = p;
b_aligned = zeros(1,cnt_parts);
b_aligned(p) = 1;
pre_v = [1,0,0];
% figure();
for i = 2:cnt_parts
    pxyz = xyz{p};
%     scatter3(pxyz(:,1),pxyz(:,2),pxyz(:,3), 15, 'filled'); hold on;
    mean_xyz = mean(pxyz);

    D = zeros(1,cnt_parts);
    min_dist = -1;
    min_idx = -1;
    % ���㵽ȫ��part�ľ���
    % ���ڵ�ǰ�ȣ�����˳ʱ�뷽����ھ�
    for j = 1:cnt_parts
        xyz2 = xyz{j};
        mean_xyz2  = mean(xyz2);
        dir = mean_xyz2-mean_xyz;
        dir = dir./norm(dir);
        if i <=2 && dot(dir(:,[1,3]),pre_v(:,[1,3]))<0 %���򲻶�
            D(j) = -1;
            continue;
        end
        if b_aligned(j) > 0%�ҵ��Ѿ����ҹ�����
            D(j) = -1;
            continue;
        end
        % ����leg�ľ���
        [D_idx,dist] = knnsearch(xyz2(end/2:end,:),pxyz(end/2:end,:),'K', 1);
%         xyz2 = xyz2(D_idx,:);
%         delta = xyz2-pxyz;
        D(j) = sum(dist);
        if min_dist == -1 || (min_dist > D(j)&&D(j)>0 )
            min_dist = D(j);
            min_idx = j;
        end
    end
    pre_v = mean(xyz{min_idx})-mean(pxyz);
    pre_v = pre_v./norm(pre_v);
    p = min_idx;
    aligned_idx(i) = p;
    b_aligned(p) = 1;
end
if aligned_idx(end) == 0
    i = find(b_aligned == 0)
    aligned_idx(end) = i;
end
aligned_xyz = xyz(aligned_idx);


