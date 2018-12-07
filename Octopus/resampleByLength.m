function [new_pts, len] = resampleByLength(xyz, n)
%% xyz: 3D coordinate
%  length: 每个点到最开始的点的距离， 第一个数是第一个点到第2个点的距离
%  n: 采样多少个点，包括开头和结尾
%% test
tmp = xyz(2:end,:)-xyz(1:end-1,:);
lengths = sqrt(dot(tmp,tmp, 2));
lengths = cumsum(lengths)';
len = lengths(end);
avg_len = max(lengths);

plength = avg_len/(n-1);
tmp_lens = 0;%会随着迭代改变
l = 1;%记录每次的关键点
new_pts = zeros(n-1, 3);
new_pts(1, :) = xyz(1,:);
%开始重新计算点的位置
for i = 2:n
    tmp_lens = (i-1)*plength;%当前要达到的长度
    %第一个大于要求length的区段， 即目标点所在的区段
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
    v = delta * v./norm(v); %应该取两个点中的哪个位置
    
    pnew = xyz(idx1,:)+v;
    new_pts(i,:) = pnew;  
end

% figure();
% plot3(new_pts(:,1), new_pts(:,2), new_pts(:,3), 'r.-'); hold on;
% plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'b.-');
