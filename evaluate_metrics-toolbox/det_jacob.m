function D=det_jacob(ddf,lab)
% DETJAC Calculates the determinants of the Jacobian matrices of a
% deformation field.
%   D = DETJAC(J) calculates the determinants of the Jacobian matrices of
%   the deformation field J, where J is a matrix with size M-by-N-by-2-by-2.
%   The output D is a matrix with size M-by-N.
% lab_sz=size(lab,4);
D=zeros(size(lab));
for lab_id=1:size(lab,4)
    J=cal_jac_vol(ddf.*lab(:,:,:,lab_id));
    % Calculate the determinants of the Jacobian matrices
    a = J(:,:,:,1,1);
    b = J(:,:,:,1,2);
    c = J(:,:,:,1,3);
    d = J(:,:,:,2,1);
    e = J(:,:,:,2,2);
    f = J(:,:,:,2,3);
    g = J(:,:,:,3,1);
    h = J(:,:,:,3,2);
    i = J(:,:,:,3,3);
    D(:,:,:,lab_id) = (a.*e.*i)+(b.*f.*g)+(c.*d.*h)-(g.*e.*c)-(h.*f.*a)-(i.*d.*b); 
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


function jacob=cal_jac_vol(vol)
vol_sz=size(vol);
% jacob=zeros([vol_sz,vol_sz(end)])
jacob=[];
for v=1:vol_sz(end)
[Gx, Gy, Gz] = imgradientxyz(vol(:,:,:,v));
G=cat(4,Gx,Gy,Gz);
jacob=cat(5,jacob,G);
end
% jacob=cat(5,Gx,Gy,Gz);
for i = 1:vol_sz(end)
jacob(:,:,:,i,i)=jacob(:,:,:,i,i)+1;
end
end

% 
% def eval_detJ_lab(vol1=None,vol2=None,disp=None,thresh=0.5):
%     ndims=3
%     label=vol1>thresh
%     label=label*(ndimage.laplace(label) < 0.1)
%     # label=vol2[...,0]>thresh
%     # label=np.ones_like(vol2[...,0])
%     # a=np.stack(np.gradient(disp,axis=[-2,-3,-4]),-1)
%     # b=np.sum(label)
%     # rescale_factor=2
%     rescale_factor=1
%     # label=zoom(label, [1./rescale_factor]*ndims+[1], mode='nearest')
%     label=label[...,::rescale_factor,::rescale_factor,::rescale_factor]
%     # disp=zoom(disp, [rescale_factor]*ndims+[1], mode='nearest')
%     # label=zoom(label, rescale_factor, mode='nearest')
%     # a=np.stack(np.gradient(disp,axis=[-2,-3,-4]),-1)
%     # b = np.linalg.det(a)
%     Jacob=np.stack(np.gradient(disp,axis=[-4,-3,-2]),-1)
%     Jacob[..., 0, 0] = Jacob[..., 0, 0] + 1
%     Jacob[..., 1, 1] = Jacob[..., 1, 1] + 1
%     Jacob[..., 2, 2] = Jacob[..., 2, 2] + 1
%     return np.sum((np.linalg.det(Jacob)<0)*label)