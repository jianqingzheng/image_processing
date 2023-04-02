function asd=avg_surf_dist(vol1,vol2,side_idx,sample_rate)
if nargin<4
    sample_rate=1;
if nargin<3
    side_idx=1;
end
end
% % lap_k1=[];
% lap_k0=[0,0,0;
%         0,-1,0;
%         0,0,0];
% lap_k2=[0,-1,0;
%         -1,6,-1;
%         0,-1,0];
% lap_kernel=cat(3,lap_k0,lap_k2,lap_k0);
% v1=convn(double(vol1),lap_kernel,'same')>0.5;
% v2=convn(double(vol2),lap_kernel,'same')>0.5;
v1=lap_conv3d(vol1)>0.5;
v2=lap_conv3d(vol2)>0.5;
% volshow(v2);
[~,coord1]=img2ind(v1);
[~,coord2]=img2ind(v2);
dist_mat=pdist2(coord1(1:sample_rate:end,:),coord2(1:sample_rate:end,:),"euclidean");
asd=mean(min(dist_mat,[],side_idx));

end

function [ind,pnd]=img2ind(img,thresh)
if nargin<2
    thresh=0.000000001;
end
img_sz=size(img);
ind=find(img>thresh);
pnd=[];
id=ind;
for s=img_sz
    coord=mod(id,s);
    pnd=[pnd,coord];
    id=floor(id/s);
end
end