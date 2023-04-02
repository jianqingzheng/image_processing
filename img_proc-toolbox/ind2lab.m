function lab=ind2lab(ind,lab_sz)
ind_tmp=round(ind);
ind_tmp(any(ind_tmp<1,2),:)=[];
ind_tmp(any(ind_tmp>lab_sz(1),2),:)=[];
lab=zeros(lab_sz);%
if isempty(ind_tmp)
    disp('is empty');
else
    lab(ind_tmp(:,1)+(ind_tmp(:,2)-1).*lab_sz(1))=1;
end
end



