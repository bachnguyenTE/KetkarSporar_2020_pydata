%%% 60s FFF analysis

addpath(genpath('your_path'))

clc;
clear all;
close all;

% CONSTANTS
REGION = 'AT';  % which group or layer of ROIs will be included in this plot
COLOR_OF_PLOT = [0 .5 0];% GREEN
COLOR_CLOUD =[.5 .7 .5];% GREEN
gen_str = GENOTYPE_STRING;

% % % % %
% pdatapath='/Users/ksporar/Documents/MATLAB/Katja_pData';
pdatapath='/your_dataPath';

addpath(pdatapath);
cd(pdatapath);

database_select_samples_ks;

% L3
indsSummaryFile  = find(i_fff60s.*i_UASGCaMP6F_L3MH56Gal4.*~i_moving); %paper


% create structures by neuron
neurStructs = create_neuron_structure_all_ks(indsSummaryFile);
neurData = load_neuron_data10Hz_byRegion(neurStructs,pdatapath, REGION);
data = aggregate_fff60sonly_means10Hz_triggeronON(neurData);


iRATE = 10; %rate at which data are interpolated
dur = size(data.rats,2);
DURS=60;

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


cur_t = [1:size(cur_mat,2)]/iRATE;


Q = corr(mean(data.stims)',cur_mat');
inds = find(Q>0.5); %finds the cells whose response correlates well with the stimulus 
cur_mat_pos = cur_mat(inds,:);
cur_IDs_pos = cur_IDs(inds);

%calculating means across flies
[x_pos,m_pos,e_pos] = mean_cat_full(cur_mat_pos,1,cur_IDs_pos);
%calculating the max for each trace
max_pos = max(x_pos,[],2)-1;


%negative correlation with stimulus
Q = corr(mean(data.stims)',cur_mat');
inds = find(Q<-0.5);
size(inds)
cur_mat_neg = cur_mat(inds,:);
cur_IDs_neg = cur_IDs(inds);

%calculating means across flies
[x_neg,m_neg,e_neg] = mean_cat_full(cur_mat_neg,1,cur_IDs_neg);
max_neg = max(x_neg,[],2)-1;


% plot negatively correlated cells, mean across flies
figure; hold on
t = [1:dur]/10;
subplot(2,1,1);
cm=colormap('lines');
h2 = plot_err_patch_v2(cur_t,m_neg-1,e_neg,[0 0.5 0],[0 0.80 0]); %m_neg-1 is there to put the calcium trace down to zero
title([gen_str ', neg corr, mean across flies']);
legend([h2],sprintf([gen_str ', N = %d ( %d )'],size(x_neg,1),size(cur_mat_neg,1)),...
    'location','northeast');
plot(cur_t, (round(mean(data.stims))*0.1)+0.6)
ylim([-0.7 1.0]);
xlim([0 140])
line([0 4],[0 0],'color',[0 0 0]); %for 2 s


subplot(2,1,2);
plot(t,x_neg)
title('individual fly means');
ylim([-0.8 3.0]); %for L2
line([DURS DURS],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 DURS],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;

%% plot negatively correlated cells, mean across ROIs
mROI_neg = mean(cur_mat_neg);
eROI_neg = std(cur_mat_neg,[],1)/sqrt(size(cur_mat_neg,1)); %S.E.M

figure; hold on
subplot(2,1,1);
cm=colormap('lines');
h2 = plot_err_patch_v2(cur_t,mROI_neg-1,eROI_neg,[0.75 0 0.75],[0.75 0.5 0.75]);
title([gen_str ', neg corr, mean by ROI']);
legend([h2],sprintf([gen_str ', N = %d ( %d )'],size(x_neg,1),size(cur_mat_neg,1)),...
    'location','northeast');
plot(cur_t, (round(mean(data.stims))*0.1)+0.6)
ylim([-0.8 2.0]); %for L2
% ylim([-0.8 1.0]); %for L3
% ylim([-0.5 0.5]); %DE
% ylim([-0.5 0.4]); %ort rescue
% line([0 10],[0 0],'color',[0 0 0]); %for 5 s
line([0 4],[0 0],'color',[0 0 0]); %for 2 s

subplot(2,1,2);
plot(t,cur_mat_neg)
title('individual fly means');
ylim([-0.8 3.0]); %for L2
% ylim([-0.8 1.0]); %for L3
% ylim([-0.5 0.5]); %DE
% ylim([-0.5 0.4]); %ort rescue
line([DURS DURS],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 DURS],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;


