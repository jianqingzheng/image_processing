function out_vol=keep_remove_ranked_vol(vol,rank_nums)
out_vol=vol>0.5;
if ndims(vol)>2
    L = bwconncomp(vol,26);% 
else
    L = bwconncomp(vol,8);%
end
stats = regionprops(L,'Area');
Ar = cat(1, stats.Area);
[~,ind]=sort(Ar,'descend');
LM = labelmatrix(L);
%     out_vol(find(LM==ind(remove_vol_rank)))=0;%
if rank_nums(1)>0
    out_vol(~ismember(LM,ind(rank_nums)))=0;%
else
    rank_nums=-1.*rank_nums;
    out_vol(ismember(LM,ind(rank_nums)))=0;%
end
end

