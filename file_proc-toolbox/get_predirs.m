function [predirs,foldernames]=get_predirs(dirs_in,layernumber)
predirs=dirs_in;
predirs_temp=predirs;
for i=0:layernumber
    if i==layernumber
        foldernames_temp=cell(size(dirs_in));
    end
    for j=1:numel(predirs)
        [predir,foldername,~]=fileparts(predirs_temp{j});
        if i==layernumber
            foldernames_temp{j}=foldername;
        else
            predirs_temp{j}=predir;
        end
    end    
end
if layernumber>0
    predirs=predirs_temp;
    
end
foldernames=foldernames_temp;
end