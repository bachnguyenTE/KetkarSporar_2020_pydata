%%% FFF with different contrast steps and gray background

addpath(genpath('/your_path'))

clc;
clear all;
close all;

% CONSTANTS

REGION = 'AT';  % which group or layer of ROIs will be included in this plot
GENOTYPE_STRING = ['L2 >> GCaMP6f, region: ' REGION];
GENOTYPE_STRING = ['L3 >> GCaMP6f, region: ' REGION];
%         GENOTYPE_STRING = ['lexAopGCaMP6fp10attp5 16HO3lexA,lexAopmyrtdTom 11-164Gal4 ermHA, region: ' REGION];
%     GENOTYPE_STRING = ['lexAopGCaMP6fp10attp5 16HO3lexA,lexAopmyrtdTom 11-164Gal4, region: ' REGION];
COLOR_OF_PLOT = [0 .5 0];% GREEN
COLOR_CLOUD =[.5 .7 .5];% GREEN
gen_str = GENOTYPE_STRING;


% % % % % 
pdatapath='/your_dataPath';
% pdatapath='/Users/ksporar/manuscript_code_figures/Figure1/Figure1D/L2_pData';

addpath(pdatapath);
cd(pdatapath);

database_select_samples_ks;

%L3
% indsSummaryFile  = find(i_fff_01steps_5sONgrayOFF.*i_UASGCaMP6F_L3MH56Gal4.*~i_moving);

%L2
indsSummaryFile  = find(i_fff_01steps_5sONgrayOFF.*i_UASGCaMP6F_L221DhhGal4.*~i_moving);


% % % % % % % % 
% create structures by neuron
neurStructs = create_neuron_structure_all_ks(indsSummaryFile);

neurData = load_neuron_data10Hz_byRegion(neurStructs,pdatapath, REGION);

data = aggregate_standingStripe_ks_grey_baseline(neurData);

iRATE = 10; %rate at which data are interpolated
dur = size(data.rats,3);

% start with all the data, keep the ones that are positively correlated with stimulus (by value Q)
cur_mat = data.rats;
cur_IDs = data.flyID;

cur_t = [1:size(cur_mat,3)]/iRATE;

counter = 0; % starting counter value!

for ii = 1:10 
    epochRats = squeeze(cur_mat(:,ii,:));
    [x,m,e] = mean_cat_full(epochRats,1,cur_IDs); %m=mean; e=std
    
    counter = counter + 1;
    
    if ismember(ii,1:5) %negative contrasts
    stimulus = data.stimstruct;
    elseif ismember(ii,6:10) %positive contrasts
    stimulus = data.stimstruct.*(-1);
    end
    
    %positive correlation with stimulus    
    Q = corr(stimulus,epochRats')
    inds = find(Q>0.1);
    cur_mat_pos = epochRats(inds,:)
    cur_IDs_pos = cur_IDs(inds)
    %calculating means across flies
    [x_pos,m_pos,e_pos] = mean_cat_full(cur_mat_pos,1,cur_IDs_pos);
    %calculating the max for each trace
    max_pos = max(x_pos,[],2)-1;


    %negative correlation with stimulus
    Q = corr(stimulus,epochRats')
    inds = find(Q<-0.1);
    size(inds)
    cur_mat_neg = epochRats(inds,:)
    cur_IDs_neg = cur_IDs(inds)
    %calculating means across flies
    [x_neg,m_neg,e_neg] = mean_cat_full(cur_mat_neg,1,cur_IDs_neg);
    max_neg = max(x_neg,[],2)-1;

    
    figure(ii); hold on
    cm=colormap('lines');
    h1 = plot_err_patch_v2(cur_t,m,e,[0 0 1],[0.5 0.5 1]);
% % %     title('some genotype, all cells');
    if (ii==1)
        title('-100% OFFgrey')
    elseif (ii==2)
        title('-80% OFFgrey')
    elseif(ii==3)
        title('-60% OFFgrey')
    elseif(ii==4)
        title('-40% OFFgrey')
    elseif(ii==5)
        title('-20% OFFgrey') 
    elseif(ii==6)
        title('20% ONgrey')
    elseif(ii==7)
        title('40% ONgrey')
    elseif(ii==8)
        title('60% ONgrey')
    elseif(ii==9)
        title('80% ONgrey')
    elseif(ii==10)
        title('100% ONgrey') 
    end
    
    if ismember(ii,1:5) %the same if you would write esleif
    peak_max(counter) = max(m(90:200)); %because for -40% it finds max in the first 1s
    elseif ismember(ii,6:10)
    peak_max(counter) = min(m);
    peak_max(counter) = max(m(140:200));
    end

    legend([h1],sprintf('N = %d ( %d )',size(x,1),size(epochRats,1)),...
        'location','northeast');
    plot(cur_t, data.stimstruct*0.1+2.25, 'k')
    ylim([-0.8 2.0]);
    line([0 5.5],[0 0],'color',[0 0 0]);

end

contrast = [-100 -80 -60 -40 -20 20 40 60 80 100]

figure
plot(contrast, peak_max,'k','Linewidth',4)
box off
ylim([-0.5 2.5]);



%% Plotting - per fly

THRESHOLD = 0;
fprintf('THRESHOLD for Signal sum: %d \n',THRESHOLD );

close all
THRESHOLD = 0;
fprintf('THRESHOLD for Signal sum: %d \n',THRESHOLD );

% take a mean of all assigned values, ignore 0
crossing_rats = data.rats;
crossing_rats(crossing_rats==0) = NaN;

% Mean over all ROI's
mean_val = squeeze(mean(crossing_rats,1,'omitnan'));


num_steps = size(mean_val,1);
numROIs = size(crossing_rats,1);

mean_calcium = zeros(1,num_steps);
err_calcium = zeros(1,num_steps);


% When are we going from grey to OFF?
on_step = [0 0 0 0 0 1 1 1 1 1];
on_step = logical(on_step);
off_step = ~ on_step;


% Mean over flies
cur_mat = data.rats;
cur_IDs = data.flyID;

%calculating means across flies
[x_neg,m_neg,e_neg] = mean_cat_full(cur_mat_neg,1,cur_IDs_neg);


% Calculate the height of the off peak

for kk =1:num_steps
  
        
    traces = squeeze(data.rats(:,kk,:)); %because mean_cat_full works only with 2 dimensions %taking all ROIs, the epoch we are on, full lenth (200)
    [x,m,e] = mean_cat_full(traces,1,data.flyID); %m=mean; e=std
    
    
    %ks: defining the time where we look for the baseline
    baseline = x(:,round(0.3*end):round(0.4*end)); %7.5s till 10s
    baselineONtogrey = x(:,round(0.5*end):round(0.6*end)); %13.75s till 15s


    %ks: defining the time where we look for the response
    search_peak = x(:,round(0.4*end)+1:round(0.8*end));  %10s till 15s
    search_peakONtogrey = x(:,round(0.5*end)+1:round(0.9*end));  %12.5s till 22.5s
    calcium = squeeze(nanmean(search_peak,1));
    calciumONtogrey = squeeze(nanmean(search_peakONtogrey,1));
    

    if on_step(kk)
        baseline = baselineONtogrey;
        peak = max(calciumONtogrey);
    else
        baseline = baseline;
        peak = max(calcium);
    end


    not_NaN_flies = size(x,1); % Fabi: would crash otherwise, we have only 8 flies (not 12) evaluated here !
    
    calciumStep = zeros(not_NaN_flies,1);

    
    calciumStep = zeros(not_NaN_flies,1);
    for ii = 1:not_NaN_flies
        if on_step(kk)
            calciumStep(ii) = max(search_peakONtogrey(ii,:)) - max(baselineONtogrey(ii,:));
        else
            calciumStep(ii) = max(search_peak(ii,:)) - max(baseline(ii,:));
        end

    end
        
    
    all_steps_all_means_flies(:,kk) = calciumStep; %use this for statistics!
    
    mean_calcium(kk) = mean(calciumStep);

    err_calcium(kk) = std(calciumStep,'omitnan')/sqrt(length(calciumStep));

    
  end  


% Plot ALL steps-------------------------------------------------
bar_all = figure;
grid on;

% Plot the data only for OFF steps
y1 = mean_calcium(off_step);
x1 = find(off_step);
e1 = err_calcium(off_step);

hold on;
% bar(x1,y1,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.0);
bar(x1,y1,'FaceColor',[0 1 1]);
errorbar(x1,y1,e1,'rx');

% Plot the data only for ON steps
y2 = mean_calcium(on_step); %ks
x2 = find(on_step);
e2 = err_calcium(on_step);

% bar(x2,y2,'FaceColor',[.5 .5 0],'EdgeColor',[.9 .9 0],'LineWidth',1.0);
bar(x2,y2,'FaceColor',[1 0 1]);
legend({'grey --> OFF','Error','ON --> grey'},'Location','northwest');
errorbar(x2,y2,e2,'rx');
% ylim([0 1.5]); %abs 
ylim([-0.5 2.0]) %no abs

title('max response');

% Add titles
ax = gca;
ax.XTick = 1:length(on_step);
ax.XTickLabels = {'grey --> 100% OFF';'grey --> 80% OFF'; 'grey --> 60% OFF';'grey --> 40% OFF';'grey --> 20% OFF';...
    '20% ON --> grey';'40% ON --> grey';'60% ON --> grey' ;'80% ON --> grey'; ...
    '100% ON --> grey'};
ax.XTickLabelRotation = 45;
ylabel('Absolute height max (OFF) peak') ;

hold off;



%% Statistics NEW!!!

for kk =1:num_steps
  
        
    traces = squeeze(data.rats(:,kk,:)); %because mean_cat_full works only with 2 dimensions %taking all ROIs, the epoch we are on, full lenth (200)
    [x,m,e] = mean_cat_full(traces,1,data.flyID); %m=mean; e=std
    
    
        %ks: defining the time where we look for the baseline
    baseline = x(:,round(0.3*end):round(0.4*end)); %7.5s till 10s
    baselineONtogrey = x(:,round(0.5*end):round(0.6*end)); %13.75s till 15s
%     baseline = squeeze(baseline);

    %ks: defining the time where we look for the response
    search_peak = x(:,round(0.4*end)+1:round(0.8*end));  %10s till 15s
    search_peakONtogrey = x(:,round(0.5*end)+1:round(0.9*end));  %12.5s till 22.5s
%     response = squeeze(response);
    calcium = squeeze(nanmean(search_peak,1));
    calciumONtogrey = squeeze(nanmean(search_peakONtogrey,1));
    
    
    

    if on_step(kk)
        baseline = baselineONtogrey;
        peak = max(calciumONtogrey);
    else
        baseline = baseline;
        peak = max(calcium);
    end


    not_NaN_flies = size(x,1); % Fabi: would crash otherwise, we have only 8 flies (not 12) evaluated here !
    
    
    calciumStep = zeros(not_NaN_flies,1);
    for ii = 1:not_NaN_flies
        if on_step(kk)
            calciumStep(ii) = max(search_peakONtogrey(ii,:)) - max(baselineONtogrey(ii,:));
        else
            calciumStep(ii) = max(search_peak(ii,:)) - max(baseline(ii,:));
        end

    end
        

    
    all_steps_all_means_flies(:,kk) = calciumStep; %use this for statistics!
    
    mean_calcium(kk) = mean(calciumStep);
    
    err_calcium(kk) = std(calciumStep,'omitnan')/sqrt(length(calciumStep));
%     err_calcium(kk) = std(calcium,'omitnan')/sqrt(length(calcium)); %is this better than baseline???
    
  end  

  
step1 = all_steps_all_means_flies(:,1);
step2 = all_steps_all_means_flies(:,2);
step3 = all_steps_all_means_flies(:,3);
step4 = all_steps_all_means_flies(:,4);
step5 = all_steps_all_means_flies(:,5);
step6 = all_steps_all_means_flies(:,6);
step7 = all_steps_all_means_flies(:,7);
step8 = all_steps_all_means_flies(:,8);
step9 = all_steps_all_means_flies(:,9);
step10 = all_steps_all_means_flies(:,10);


%grey to 100OFF vs 100ON to grey
[h,p,ci,stats] = ttest2(step1,step10);
p_step1 = p
h_step1 = h
stats_step1 = stats

% [a p1] = ttest2(step1 , step10{2})

%grey to 80OFF vs 80ON to grey
[h,p,ci,stats] = ttest2(step2,step9);
p_step2 = p
h_step2 = h
stats_step2 = stats

%grey to 60OFF vs 60ON to grey
[h,p,ci,stats] = ttest2(step3,step8);
p_step3 = p
h_step3 = h
stats_step3 = stats

%grey to 40OFF vs 40ON to grey
[h,p,ci,stats] = ttest2(step4,step7);
p_step4 = p
h_step4 = h
stats_step4 = stats

%grey to 20OFF vs 20ON to grey
[h,p,ci,stats] = ttest2(step5,step6);
p_step5 = p
h_step5 = h
stats_step5 = stats

% L3: 1 and 2 are stat different!!! p1 = 0.0482, p2 = 0.0333, p3 = 0.066, p4 = 0.066, p5 = 0.07
% L2: none of them! p1 = 0.173, p2 = 0.194, p3 = 0.269, p4 = 0.1775, p5 = 0.417

