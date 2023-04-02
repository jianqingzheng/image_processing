function out_vol=lap_conv3d(vol)
% lap_k1=[];
lap_k0=[0,0,0;
        0,-1,0;
        0,0,0];
lap_k2=[0,-1,0;
        -1,6,-1;
        0,-1,0];
lap_kernel=cat(3,lap_k0,lap_k2,lap_k0);
out_vol=convn(double(vol),lap_kernel,'same');
end
