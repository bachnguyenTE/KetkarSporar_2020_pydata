function out=create_neuron_structure_all_ks(inds)

T=readtable('Summary.xlsx','ReadVariableNames',1);
fname=T.fname;

out=[];

count=1;
for ii=1:length(inds)
    if exist(fname{inds(ii)})
        x=load(fname{inds(ii)});
        ns=[1:size(x.strct.dRatio,1)];
        for jj=1:length(ns)
            out(count).name=fname{inds(ii)};
            out(count).cell=ns(jj);
            count=count+1;
        end
    end
end



