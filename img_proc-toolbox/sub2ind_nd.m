function ind = sub2ind_nd( tensor_sz, subs )
% sub2ind_nd - Get the index number using the given subscript coordinates
% and the tensor size
%--------------------------------------------------------------------------
%   [ind] = sub2ind_nd( tensor_sz, subs )
%   'ind'       - the output index numbers
%   'tensor_sz' - the size of the required tensor
%   'sub'       - the given subscript corrdinates
%--------------------------------------------------------------------------
%   Examples:
%       >> ind = sub2ind_nd([5,5],[3,5;2,3])
% 
%       ind =
%           23
%           12
% 
%       >> ind = sub2ind_nd([3,3,3],[1,3,1;1,2,2])
% 
%       ind =
%           7
%           13
%--------------------------------------------------------------------------
%   MATLAB Ver R2019a
%--------------------------------------------------------------------------
%   $ Author: Jachin (Jianqing Zheng) $
%   $ Revision: 1.0 $  $ Date: 2023/04/05$
%--------------------------------------------------------------------------
%   See also: 'sub2lab', 'ind2sub'

%--------------------------------------------------------------------------
subs=round(subs);
ndims=size(subs,2);
ind=1;
tensor_sz=[1,tensor_sz];
for i=1:ndims
    ind=ind+(subs(:,i)-1).*prod(tensor_sz(1:i));
end
end