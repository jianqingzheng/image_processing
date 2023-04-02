function cd=centroid_dist(vol1,vol2)
s = regionprops3(vol1,"Centroid","PrincipalAxisLength");
coord1 = s.Centroid;
s = regionprops3(vol2,"Centroid","PrincipalAxisLength");
coord2 = s.Centroid;
% v1=lap_conv3d(vol1)>0.5;
% v2=lap_conv3d(vol2)>0.5;
% volshow(v2);
% [~,coord1]=img2ind(v1);
% [~,coord2]=img2ind(v2);
dist_mat=pdist2(coord1,coord2,"euclidean");
cd=mean(min(dist_mat,[],1));

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