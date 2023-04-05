function lab=sub2lab(sub,lab_sz)
% sub2lab - Get a binary volume (label) with one values at the subcript
% position and zero at the others
%--------------------------------------------------------------------------
%   [lab] = sub2lab(sub,lab_sz)
%   'lab'   - output binary volume with one values at the required locations
%   'sub'   - the given subscript corrdinates
%   'lab_sz'- the size of the required volume
%--------------------------------------------------------------------------
%   Examples:
%       >> lab = sub2lab([1,1;2,2],[3,3])
% 
%       lab =
%           1     0     0
%           0     1     0
%           0     0     0
% 
%       >> lab = sub2lab([1,1,1;2,2,2],[3,3,3])
% 
%       lab(:,:,1) =
%           1     0     0
%           0     0     0
%           0     0     0
% 
%       lab(:,:,2) =
%           0     0     0
%           0     1     0
%           0     0     0
% 
%       lab(:,:,3) =
%           0     0     0
%           0     0     0
%           0     0     0
%--------------------------------------------------------------------------
%   MATLAB Ver R2019a
%--------------------------------------------------------------------------
%   $ Author: Jachin (Jianqing Zheng) $
%   $ Revision: 1.0 $  $ Date: 2023/04/05$
%--------------------------------------------------------------------------
%   See also: 'ind2lab', 'lab2ind', 'sub2ind_nd', 'ind2sub'

%--------------------------------------------------------------------------
sub_tmp=round(sub);
for d=1:numel(lab_sz)
sub_tmp(any(sub_tmp(d)<1,2),:)=[];
sub_tmp(any(sub_tmp>lab_sz(d),2),:)=[];
end
lab=zeros(lab_sz);%
if isempty(sub_tmp)
    disp('is empty');
else
    lab(sub2ind_nd(lab_sz,sub_tmp))=1;
end
end



