function lab=ind2lab(ind,lab_sz)
% ind2lab - Get a binary volume (label) with one values at the indexed
% position and zero at the others
%--------------------------------------------------------------------------
%   [lab] = ind2lab(ind,lab_sz)
%   'lab'   - output binary volume with one values at the required locations
%   'ind'   - the given indexed location
%   'lab_sz'- the size of the required volume
%--------------------------------------------------------------------------
%   Examples:
%      >> lab = ind2lab([2,6],[3,3])
% 
%      lab =
%           0     0     0
%           1     0     0
%           0     1     0
%--------------------------------------------------------------------------
%   MATLAB Ver R2019a
%--------------------------------------------------------------------------
%   $ Author: Jachin (Jianqing Zheng) $
%   $ Revision: 1.0 $  $ Date: 2023/04/05$
%--------------------------------------------------------------------------
%   See also: 'sub2lab', 'lab2ind', 'sub2ind_nd'

%--------------------------------------------------------------------------
ind_tmp=round(ind);
ind_tmp(any(ind_tmp<1,2),:)=[];
ind_tmp(ind_tmp>prod(lab_sz))=[];
lab=zeros(lab_sz);%
if isempty(ind_tmp)
    disp('is empty');
else
    lab(ind_tmp)=1;
end
end



