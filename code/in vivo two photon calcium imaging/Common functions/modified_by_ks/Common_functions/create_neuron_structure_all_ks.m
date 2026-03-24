function out=create_neuron_structure_all_ks(inds)

T=readtable('Summary.xlsx','ReadVariableNames',1);
fname=T.fname;

out=[];

count=1;
for ii=1:length(inds)
    %     ns=str2num(activecells{inds(ii)});
    if exist(fname{inds(ii)})
        %         if (length(ns)==0)
        x=load(fname{inds(ii)});
        ns=[1:size(x.strct.dRatio,1)];
        %         end
        for jj=1:length(ns)
            out(count).name=fname{inds(ii)};
            out(count).cell=ns(jj);
            count=count+1;
        end
    end
end



