function out = aggregate_fff60sonly_means10Hz_triggeronON(in)

dxi_all = zeros(length(in),1);
%determine the epoch length
for ii=1:length(in)
    ii
        s=in(ii).istim;
        xi=[1 find(diff(s(1:end))<0)+1] % trigger on going to dark...
        dxi = mean(diff(xi));
        dxi_all(ii,:)=dxi;
end
epochlength = round(mean(dxi_all));
% dur = 1.3*epochlength; %works only for L2
dur = epochlength; 



rats=zeros(length(in),dur);
stims=rats;

mr=zeros(length(in),dur);
ms=zeros(length(in),dur);
name = cell(length(in),1);
for ii=1:length(in)
        s=in(ii).istim;
        xi=[1 find(diff(s(1:end-dur))<0)+1];
        iratio = in(ii).iratio;
        iratio = iratio./mean(iratio) - 1;
        d_iratio = iratio+1;
    for jj=1:length(xi)
        temp = d_iratio(xi(jj):xi(jj)+dur-1);%/median(d_iratio(xi(jj):xi(jj)+99));
        mr(jj,:)=temp;
        temp=in(ii).istim(xi(jj):xi(jj)+dur-1);
        temp=temp-min(temp); temp=temp/(max(temp));
        ms(jj,:)=temp;
    end

    rats(ii,:)=mean(mr(1:length(xi),:),1); %averaging over stimulus repititions
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
numSamplesPerMean=epochlength/2;
out.index=mean(rats(:,1:epochlength/2),2)-mean(rats(:,epochlength/2+1:epochlength),2);
ratvar=(var(rats(:,1:epochlength/2),[],2)/numSamplesPerMean+var(rats(:,epochlength/2+1:epochlength),[],2)/numSamplesPerMean).^.5; % approx 50 samples in each mean...
out.index_z=out.index./ratvar;

% Identify neurons with an inverse response
% Rdiff=mean(rats(:,1:200),2)-mean(rats(:,201:400),2);
f_norm=find(out.index_z>0.5);
f_inv=find(out.index_z<-0.5);

ufly = unique(flyID);
meanbyfly_inv = zeros(length(ufly),dur);
meanbyfly_norm = zeros(length(ufly),dur);
inv_inds = zeros(length(ufly),1);
norm_inds = zeros(length(ufly),1);
for ii=1:length(ufly)
    f=find((flyID==ufly(ii)));
    ff=intersect(f_norm,f);
    if(~isempty(ff))
        meanbyfly_norm(ii,:)=mean(rats(ff,:),1);
        norm_inds(ii) = 1;
    else
        disp('empty');
    end
    ff=intersect(f_inv,f);
    if(~isempty(ff))
        meanbyfly_inv(ii,:)=mean(rats(ff,:),1);
        inv_inds(ii) = 1;
    end
    meanbyfly_all(ii,:)=mean(rats(f,:),1); %ms: unbiased, all cells
end

meanbyfly_norm = meanbyfly_norm(norm_inds==1,:);
meanbyfly_inv = meanbyfly_inv(inv_inds==1,:);

allm_norm=mean(meanbyfly_norm,1);
alle_norm=std(meanbyfly_norm,[],1)/sqrt(sum(norm_inds));

allm_inv=mean(meanbyfly_inv,1);
alle_inv=std(meanbyfly_inv,[],1)/sqrt(sum(inv_inds));

allm_all=mean(meanbyfly_all,1);%ms
alle_all=std(meanbyfly_all,[],1)/sqrt(sum(norm_inds)); %ms: what does [] do here?

%means only across cells (by fly), not across flies yet
out.meanbyfly_norm=meanbyfly_norm;
out.meanbyfly_inv=meanbyfly_inv;
out.meanbyfly_all=meanbyfly_all;%ms: added
%mean and sem across flies
out.allm_norm=allm_norm;
out.allm_inv=allm_inv;
out.allm_all=allm_all; %ms: added
out.alle_norm=alle_norm;
out.alle_inv=alle_inv;
out.alle_all=alle_all; %ms: added

out.norm_rats=rats(f_norm,:);
out.inv_rats=rats(f_inv,:);
out.all_rats=rats;
