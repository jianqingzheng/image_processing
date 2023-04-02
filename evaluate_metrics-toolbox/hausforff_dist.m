function hd=hausforff_dist(vol1,vol2,side_idx)
if nargin<3
    side_idx=0;
end
v1=lap_conv3d(vol1)>0.5;
v2=lap_conv3d(vol2)>0.5;
[~,coord1]=img2ind(v1);
[~,coord2]=img2ind(v2);
if side_idx==1 || side_idx==2
dist_mat=pdist2(coord1,coord2,"euclidean");
hd=max(min(dist_mat,[],side_idx));
else
[hd]=HausdorffDist(coord1,coord2);
end
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