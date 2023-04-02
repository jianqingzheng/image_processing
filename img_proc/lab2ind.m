function [ind,pnd]=lab2ind(lab,thresh)
if nargin<2
    thresh=0.1;
end
img_sz=size(lab);
ind=find(lab>thresh);
pnd=[];
id=ind;
for s=img_sz
    coord=mod(id,s);
    pnd=[pnd,coord];
    id=floor(id/s);
end
end
