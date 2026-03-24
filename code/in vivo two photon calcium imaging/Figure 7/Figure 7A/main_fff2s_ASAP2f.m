% set up to make all fff figures...
clear all;
clc;

pdatapath='/your_datapath';
addpath(pdatapath);
cd(pdatapath);

database_select_samples_ks;

%full field flash with ASAP2f
indsSummaryFile = find(i_fff2s.*i_L221DhhGal4_Asap2f.*~i_moving.*((flyID==722)|(flyID==739)|(flyID==740)|(flyID==741)|(flyID==742)|(flyID==743)|(flyID==744)|(flyID==745)|(flyID==751)|(flyID==753)|(flyID==754)));

% indsSummaryFile = find(i_fff2s.*i_L3MH56Gal4_Asap2f.*~i_moving.*((flyID==726)|(flyID==727)|(flyID==728)|(flyID==729)|(flyID==731)|(flyID==732)|(flyID==749)|(flyID==750)|(flyID==756)|(flyID==757)|(flyID==758)|(flyID==759)|(flyID==761)|(flyID==764)|(flyID==765)));




% create structures by neuron
neurStructs = create_neuron_structure_all_ks(indsSummaryFile);

% load all the data!
% neurData = load_neuron_data10Hz_ks(fTm9LexA,pdatapath);
neurData = load_neuron_data40Hz_ks(neurStructs,pdatapath);


data = aggregate_fff300ms_means10Hz_BleedThruFix(neurData); %in = xTm9LexA


%%
%save('fff_analysis_res_Apr24','mL2_562');

% Paper figures
%save fff_analysis_L4_L2_silenced;
%cd 'C:/CalciumImaging';
%load('fff_analysis_L4_L2_silenced'); 

dur = size(data.rats,2);
durs = 0.3;

%Plotting Tm9LexA, "normal" and "inverted" cells
t = [1:dur]/10;
figure; hold on;
subplot(2,2,1);
cm=colormap('lines');
h1 = plot_err_patch_v2(t,[data.allm_norm-1 ],...
    [data.alle_norm ],[0 0 1],[0.5 0.5 1]);
h2 = plot_err_patch_v2(t ,[data.allm_inv-1 ],...
   [data.alle_inv],[0.75 0 0.75],[0.75 0.5 0.75]);
% plot([0 5],0.3*[1 1],'k-','linewidth',3);
% plot([10 15],0.3*[1 1],'k-','linewidth',3);
title('L2 ASAP2s, all "normal" and "inverted" cells');
legend([h1, h2],sprintf('ASAP2s>>GCaMP6F, "normal", N = %d ( %d )',size(data.meanbyfly_norm,1),size(data.norm_rats,1)),...
    sprintf('ASAP2s>>GCaMP6F, "inverted", N = %d ( %d )',size(data.meanbyfly_inv,1),size(data.inv_rats,1)),...
    'location','northeast');
plot(t, (round(mean(data.stims))*0.1)+0.05)

xlabel('time (sec)');ylabel('response (dR/R)');
ylim([-0.1 0.1]);
xlim([0 dur/10]);
line([0 10],[0 0],'color',[0 0 0]);
set(gca,'xTick',0:5:10);
set(gcf,'Color','w');
niceaxes;
%Plotting Tm9LexA, "normal" and "inverted" cells
t = [1:dur]/10;
figure; hold on;
subplot(2,2,1);
cm=colormap('lines');
h1 = plot_err_patch_v2(t,[data.allm_norm-1 ],...
    [data.alle_norm ],[0 0 1],[0.5 0.5 1]);
h2 = plot_err_patch_v2(t ,[data.allm_inv-1 ],...
   [data.alle_inv],[0.75 0 0.75],[0.75 0.5 0.75]);
% plot([0 5],0.3*[1 1],'k-','linewidth',3);
% plot([10 15],0.3*[1 1],'k-','linewidth',3);
title('L2 ASAP2s, all "normal" and "inverted" cells');
legend([h1, h2],sprintf('ASAP2s>>GCaMP6F, "normal", N = %d ( %d )',size(data.meanbyfly_norm,1),size(data.norm_rats,1)),...
    sprintf('ASAP2s>>GCaMP6F, "inverted", N = %d ( %d )',size(data.meanbyfly_inv,1),size(data.inv_rats,1)),...
    'location','northeast');
plot(t, (round(mean(data.stims))*0.1)+0.05)

xlabel('time (sec)');ylabel('response (dR/R)');
ylim([-0.1 0.1]);
xlim([0 dur/10]);
line([0 10],[0 0],'color',[0 0 0]);
set(gca,'xTick',0:5:10);
set(gcf,'Color','w');
niceaxes;

subplot(2,2,3);
plot(t,data.meanbyfly_norm-1)
title('normal, DAC metric');
ylim([-0.1 0.1]);
line([durs durs],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 durs+2],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;


subplot(2,2,4);
plot(t,data.meanbyfly_inv-1)
title('inverted, DAC metric');
ylim([-0.1 0.1]);
line([durs durs],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 durs+2],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;







%%
% correlate with stimulus to select normal and inverted ROIs
iRATE = 40; %rate at which data are interpolated
dur = size(data.rats,2);


% start with all the data, keep the ones that are positively correlated
% with stimulus (by value Q)
cur_mat = data.rats;
cur_IDs = data.flyID;
%clean up by zeros-only datasets and NaN dataset
inds = find(sum(cur_mat,2)~=0);
cur_mat = cur_mat(inds,:);
cur_IDs = cur_IDs(inds);
inds = ~isnan(sum(cur_mat,2));%ms: added these two lines, because there were NaN in dataset
cur_mat = cur_mat(inds,:);
cur_IDs = cur_IDs(inds);
% mean_res = mean(cur_mat,1); %was correlating with the mean before
% Q = corr(mean_res',cur_mat');
cur_t = [1:size(cur_mat,2)]/iRATE;


%positive correlation with stimulus

Q = corr(mean(data.stims)',cur_mat');
inds = find(Q>0.2); %finds the cells whose response correlates well with
% the stimulus 
cur_mat_pos = cur_mat(inds,:);
cur_IDs_pos = cur_IDs(inds);

%calculating means across flies
[x_pos,m_pos,e_pos] = mean_cat_full(cur_mat_pos,1,cur_IDs_pos);
%calculating the max for each trace
max_pos = max(x_pos,[],2)-1;
    


%negative correlation with stimulus

Q = corr(mean(data.stims)',cur_mat');
inds = find(Q<-0.2);
size(inds)
cur_mat_neg = cur_mat(inds,:);
cur_IDs_neg = cur_IDs(inds);

%calculating means across flies
[x_neg,m_neg,e_neg] = mean_cat_full(cur_mat_neg,1,cur_IDs_neg);
max_neg = max(x_neg,[],2)-1;

subplot(2,2,3);
plot(t,data.meanbyfly_norm-1)
title('normal, DAC metric');
% ylim([-0.2 0.2]); %ASAP
ylim([-0.5 0.9]);
line([durs durs],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 durs+2],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;


subplot(2,2,4);
plot(t,data.meanbyfly_inv-1)
title('inverted, DAC metric');
% ylim([-0.2 0.2]);
ylim([-0.5 0.9]);
line([durs durs],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 durs+2],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;

%%
% correlate with stimulus to select normal and inverted ROIs
iRATE = 40; %rate at which data are interpolated
dur = size(data.rats,2);


% start with all the data, keep the ones that are positively correlated
% with stimulus (by value Q)
cur_mat = data.rats;
cur_IDs = data.flyID;
%clean up by zeros-only datasets and NaN dataset
inds = find(sum(cur_mat,2)~=0);
cur_mat = cur_mat(inds,:);
cur_IDs = cur_IDs(inds);
inds = ~isnan(sum(cur_mat,2));%ms: added these two lines, because there were NaN in dataset
cur_mat = cur_mat(inds,:);
cur_IDs = cur_IDs(inds);
% mean_res = mean(cur_mat,1); %was correlating with the mean before
% Q = corr(mean_res',cur_mat');
cur_t = [1:size(cur_mat,2)]/iRATE;


%positive correlation with stimulus

Q = corr(mean(data.stims)',cur_mat');
inds = find(Q>0.2); %finds the cells whose response correlates well with
% the stimulus 
cur_mat_pos = cur_mat(inds,:);
cur_IDs_pos = cur_IDs(inds);

%calculating means across flies
[x_pos,m_pos,e_pos] = mean_cat_full(cur_mat_pos,1,cur_IDs_pos);
%calculating the max for each trace
max_pos = max(x_pos,[],2)-1;
    


%negative correlation with stimulus

Q = corr(mean(data.stims)',cur_mat');
inds = find(Q<-0.2);
size(inds)
cur_mat_neg = cur_mat(inds,:);
cur_IDs_neg = cur_IDs(inds);

%calculating means across flies
[x_neg,m_neg,e_neg] = mean_cat_full(cur_mat_neg,1,cur_IDs_neg);
max_neg = max(x_neg,[],2)-1;


% plot positively correlated cells, mean across flies
figure; hold on
cm=colormap('lines');
% h2 = plot_err_patch_v2(cur_t,m_neg-1,e_neg,[0.75 0 0.75],[0.75 0.5 0.75]);
h2 = plot_err_patch_v2(cur_t,m_pos-1,e_pos,[0 0.5 0],[0 0.80 0]);
title('L3 UASGCaMP6f');
legend([h2],sprintf('L3>>GCaMP6F, N = %d ( %d )',size(x_pos,1),size(cur_mat_pos,1)),...
    'location','northeast');
plot(cur_t, (round(mean(data.stims))*0.05)+0.03)
% ylim([-0.1 0.1]); %ASAP
ylim([-0.5 1.0]);
line([0 0.6],[0 0],'color',[0 0 0]);



% plot positively correlated cells, mean across flies
figure; hold on
subplot(2,1,1); hold on;
cm=colormap('lines');
% h2 = plot_err_patch_v2(cur_t,m_neg-1,e_neg,[0.75 0 0.75],[0.75 0.5 0.75]);
h2 = plot_err_patch_v2(cur_t,m_pos-1,e_pos,[0 0.5 0],[0 0.80 0]);
title('L3 UASGCaMP6f');
legend([h2],sprintf('L3>>GCaMP6F, N = %d ( %d )',size(x_pos,1),size(cur_mat_pos,1)),...
    'location','northeast');
plot(cur_t, (round(mean(data.stims))*0.05)+0.03)
% ylim([-0.1 0.1]); %ASAP
ylim([-0.5 1.0]);
line([0 0.6],[0 0],'color',[0 0 0]);

subplot(2,1,2);
plot(t,x_pos)
title('individual fly means');
ylim([0.9 1.1]); %for L2
% ylim([-0.8 1.0]); %for L3
% ylim([-0.5 0.5]); %DE
% ylim([-0.5 0.5]); %ort rescue
line([durs durs],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 durs],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;

%calculate mean across ROIs
mROI_pos = mean(cur_mat_pos);
eROI_pos = std(cur_mat_pos,[],1)/sqrt(size(cur_mat_pos,1)); %S.E.M

% plot by mean across ROIs
% Katja
figure; hold on
% subplot(2,2,1); hold on;
cm=colormap('lines');
h2 = plot_err_patch_v2(cur_t,mROI_pos-1,eROI_pos,[0 0.5 0],[0 0.80 0]);
title('L3 UASGCaMP6f');
legend([h2],sprintf('L3>>GCaMP6F, N = %d ( %d )',size(x_pos,1),size(cur_mat_pos,1)),...
    'location','northeast');
plot(cur_t, (round(mean(data.stims))*0.05)+0.03)
ylim([-0.1 0.1]); %ASAP
% ylim([-0.5 1.0]);
line([0 0.6],[0 0],'color',[0 0 0]);







