function out = aggregate_fffall_means10Hz_BleedThruFix_v2(in)

dxi_all = zeros(length(in),1);
%determine the epoch length
for ii=1:length(in)
        s=in(ii).istim;
        xi=find(diff(s(1:end))<0)+1; % trigger on going to dark...
        dxi = mean(diff(xi));
        dxi_all(ii,:)=dxi;
end
epochlength = round(mean(dxi_all));
dur = epochlength*1.3; %no longer *1.3

rats=zeros(length(in),dur);
stims=rats;

mr=zeros(length(in),dur);
ms=zeros(length(in),dur);
name = cell(length(in),1);
for ii=1:length(in)
    s=in(ii).istim;
        xi=find(diff(s(1:end-dur))<0)+1; % trigger on going to dark... diff(s) is -1 then
        iratio = in(ii).iratio;
        
        iratio = iratio./mean(iratio)-1; %calculating dF/F with F being the mean of the whole trace 
        d_iratio = iratio;
        
    for jj=1:length(xi)
        temp = d_iratio(xi(jj):xi(jj)+dur-1);

        mr(jj,:)=temp;
        
        temp=in(ii).istim(xi(jj):xi(jj)+dur-1);
        ms(jj,:)=temp;
    end
    rats(ii,:)=mean(mr(1:length(xi),:),1); %averaging over stimulus repetitions
    stims(ii,:)=mean(ms(1:length(xi),:),1);
    name{ii}=in(ii).name;

end
neuron = [in.cell];
flyID = [in.flyID];

out.rats=rats;
out.stims=stims;
out.neuron=neuron;
out.flyID=flyID;
out.name=name;



