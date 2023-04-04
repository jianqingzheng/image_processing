function out_vol=serial_dilate_erode_vol(vol,dilate_erode_rates,growing_msk,min_area,se_type)
if nargin<5
    se_type='sphere';
    if nargin<4
        min_area=0;
        if nargin<3
            growing_msk=1;
        end
    end
end
for r= dilate_erode_rates
    se=strel(se_type,abs(r));
%     se=strel('cube',abs(r));
    if r>0
        vol=imdilate(vol,se);
    else
        vol=imerode(vol,se);
    end
    vol=vol.*growing_msk;
    if min_area>0
        vol=bwareaopen(vol,min_area);
    end
end
out_vol=vol;
end