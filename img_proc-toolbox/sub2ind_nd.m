function ind = sub2ind_nd( vol_sz, subs )
% sz_subs=size(subs);
% ndims=numel(vol_sz);
subs=round(subs);
ndims=size(subs,2);
ind=1;
% index = cell(1, num_dims);
% index(:) = {':'};
vol_sz=[1,vol_sz];
for i=1:ndims
%     index{end} = i;
%     ind=ind+(subs(index{:})-1).*prod(vol_sz(1:i));
    ind=ind+(subs(:,i)-1).*prod(vol_sz(1:i));
end
end