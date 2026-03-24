%% 5s FFF analysis
% you will need to specify paths on line 8 and line 46, and select lines 54 or 55 for L2 or L3 respectively.

clc;
clearvars -EXCEPT slob* contr*;
close all;

addpath(genpath('/your_path'))

% CONSTANTS

REGION = 'AT';  % which group or layer of ROIs will be included in this plot
gen_str = ['L2 >> GCaMP6f, region: ' REGION];
COLOR_OF_PLOT = [0 .5 0];% GREEN
COLOR_CLOUD =[.5 .7 .5];% GREEN

% % % % % 
pdatapath='/your_datapath';

addpath(pdatapath);
cd(pdatapath);

database_select_samples_ks;


% indsSummaryFile  = find(i_fff5s.*i_UASGCaMP6F_L221DhhGal4_cross_to_w.*~i_moving);
indsSummaryFile  = find(i_fff5s.*i_UASGCaMP6f_L3MH56Gal4_cross_to_w.*~i_moving);

% create structures by neuron
neurStructs = create_neuron_structure_all_ks(indsSummaryFile);
neurData = load_neuron_data10Hz_byRegion(neurStructs,pdatapath, REGION);
data = aggregate_fffall_means10Hz_BleedThruFix_v2_manuscript(neurData);


%%
% correlate with stimulus to select normal and inverted ROIs
iRATE = 10; %rate at which data are interpolated
dur = size(data.rats,2);
DURS = 2;


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


%positive correlation with stimulus
Q = corr(nanmean(data.stims)',cur_mat');
% Q = corr(mean(mTm9LexA.stims)',cur_mat');
inds = find(Q>0.5); %finds the cells whose response correlates well with the stimulus 
cur_mat_pos = cur_mat(inds,:);
cur_IDs_pos = cur_IDs(inds);

%calculating means across flies
[x_pos,m_pos,e_pos] = mean_cat_full(cur_mat_pos,1,cur_IDs_pos);
%calculating the max for each trace
max_pos = max(x_pos,[],2)-1;


%negative correlation with stimulus
Q = corr(nanmean(data.stims)',cur_mat');
inds = find(Q<-0.5);
size(inds)
cur_mat_neg = cur_mat(inds,:);
cur_IDs_neg = cur_IDs(inds);

%calculating means across flies
[x_neg,m_neg,e_neg] = mean_cat_full(cur_mat_neg,1,cur_IDs_neg);
max_neg = max(x_neg,[],2)-1;


% plot negatively correlated cells, mean across flies
figure; hold on
subplot(2,1,1); 
cm=colormap('lines');
h2 = plot_err_patch_v2(cur_t,m_neg,e_neg,[0 0.5 0],[0 0.80 0]); %green
title([gen_str ', neg corr, mean by fly']);
legend([h2],sprintf([gen_str ', N = %d ( %d )'],size(x_neg,1),size(cur_mat_neg,1)),...
    'location','northeast');
plot(cur_t, (round(mean(data.stims))*0.1)+0.6)
ylim([-1.0 2.5]);
line([0 4],[0 0],'color',[0 0 0]); %for 2 s


% plot negatively correlated cells, individual flies
subplot(2,1,2);
t = [1:dur]/10;
plot(t,x_neg)
title('individual fly means');
ylim([-1.0 2.5]);
line([DURS DURS],[-0.6 0.6],'color',[0 0 0],'linestyle','--');
line([0 DURS],[0 0],'color',[0 0 0]);
set(gcf,'Color','w');
niceaxes;


