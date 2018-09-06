function [seg, pend] = tag(seg, dir, idx_n, ii, t)
    pend = ii;
    if seg(ii) ~= 0 
        return ; 
    end

    seg(ii) = t;
    init_nei = idx_n{ii};%列向量
    init_nei = init_nei(:,1);
    tmp = find(seg(init_nei) ~= t);
    idx_nei = init_nei(tmp);

%     %debug
%     if ~isempty([find(init_nei==10, 1),find(init_nei==15, 1)])
%         init_nei; 
%     end
%     if ii == 10
%         init_nei;
%     end
    
    if ~isempty(idx_nei)
        idx_nei = idx_nei(1); %往最近的点扩展
        dir_nei = dir(idx_nei,:);
        ang = abs(sum(dir(ii,:).*dir_nei, 2)); %列向量

        % 下一个被标记的点应该是方向满足要求的条件下距离最近的
        if ang > 0.707
           j = idx_nei(1);
           [seg, pend] = tag(seg, dir, idx_n, j, t);
        end
    end
   
end