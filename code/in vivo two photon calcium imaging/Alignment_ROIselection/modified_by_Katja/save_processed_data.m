function save_processed_data(in);


fileloc = in.fileloc;

curDir = pwd;

cd(in.fileloc);

strct=in;
if(isfield(in,'ch1'))
    strct = rmfield(strct,'ch1');
    strct = rmfield(strct,'ch1a');
    strct = rmfield(strct,'ch2');
    strct = rmfield(strct,'ch2a');
end

pathroot=in.fileloc;
f=find((pathroot == '\') + (pathroot == '/'));
f=f(end);
pathroot = pathroot((f+1):end); 

save([pathroot '_pData.mat'],'strct');

cd(curDir);