function out_vol=serial_dilate_erode_vol(vol,dilate_erode_rates,growing_msk)
if nargin<3
    growing_msk=1;
end
for r= dilate_erode_rates
    se=strel('sphere',abs(r));
%     se=strel('cube',abs(r));
    if r>0
        vol=imdilate(vol,se);
    else
        vol=imerode(vol,se);
    end
    vol=vol.*growing_msk;
end
out_vol=vol;
end