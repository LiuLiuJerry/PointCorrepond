function [new_pts, len] = resampleByLength(xyz, n)
%% xyz: 3D coordinate
%  length: ÿ���㵽�ʼ�ĵ�ľ��룬 ��һ�����ǵ�һ���㵽��2����ľ���
%  n: �������ٸ��㣬������ͷ�ͽ�β
%% test
tmp = xyz(2:end,:)-xyz(1:end-1,:);
lengths = sqrt(dot(tmp,tmp, 2));
lengths = cumsum(lengths)';
len = lengths(end);
avg_len = max(lengths);

plength = avg_len/(n-1);
tmp_lens = 0;%�����ŵ����ı�
l = 1;%��¼ÿ�εĹؼ���
new_pts = zeros(n-1, 3);
new_pts(1, :) = xyz(1,:);
%��ʼ���¼�����λ��
for i = 2:n
    tmp_lens = (i-1)*plength;%��ǰҪ�ﵽ�ĳ���
    %��һ������Ҫ��length�����Σ� ��Ŀ������ڵ�����
    while l<=size(lengths,2) &&  lengths(l) <= tmp_lens
        l = l+1; 
    end
    idx1 = l;
    idx2 = l+1;
    if idx2 > size(lengths,2)+1
       new_pts(i, :) = xyz(end, :);
       continue;
    end
    if l == 1
        delta = tmp_lens;
    else
        delta = tmp_lens-lengths(l-1);
    end
    v = xyz(idx2,:)-xyz(idx1,:);
    v = delta * v./norm(v); %Ӧ��ȡ�������е��ĸ�λ��
    
    pnew = xyz(idx1,:)+v;
    new_pts(i,:) = pnew;  
end

% figure();
% plot3(new_pts(:,1), new_pts(:,2), new_pts(:,3), 'r.-'); hold on;
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.-');
