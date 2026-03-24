%%
% This script evaluates the data from the OFF OFF stimulus from adapted
% background. The data is displayed by different plots. 
% 
%

clear all;
close all;
clc;
autosave = 0; %1 if you want to save the plots automatically to the savelocation defined above

pdatapath='/Users/ksporar/Documents/MATLAB/Katja_pData';
addpath(pdatapath);

addpath(pdatapath);
cd(pdatapath);
OFF_stimulus = 1; %1 in case its the OFF OFF stimulus 0 in case its the ON ON stimulus
ndFilter = 0.6;   %= 0 or 0.6
if ndFilter == 0
NDFilter = '_zeroND'; 
elseif ndFilter == 0.6
    NDFilter = '_zerosixND';
end 
decay_analysis = 0; %1 if you want to calculate the decay rate of the A step
% fetch Information from Masterfoldersummary
database_select_samples_ks;


%These are the measured luminance and contrast values for the stimulus
real_contrast = [-25 -25 -25 -25 -25 -25];
real_contrast1 = [-0.250000001875000 -0.357142858750000 -0.464285715625000 -0.571428572500000 -0.678571429375000 -0.785714286250000].*100;%[-0.262831359 -0.362098139 -0.476593345 -0.574168077 -0.688663283 -0.786238015];

real_luminance = [75 64 54 43 32 21];%[0.435666667 0.377 0.309333333 0.251666667 0.184 0.126333333 ];
real_luminance1 = [57 48 40 32 24 16];%[0.329 0.279666667 0.230666667 0.184 0.135666667 0.096666667 ];
  
Difference_luminance = [-0.155333333, -0.214 -0.281666667 -0.339333333 -0.407 -0.464666667];
Difference_luminance1 =[-0.106666667 -0.097333333 -0.078666667 -0.067666667 -0.048333333 -0.029666667];

contrast_baseline_bstep = [-0.443316413 -0.52679075 -0.609701072 -0.688663283... 
-0.770445572 -0.83643542];


OS = computer;
if OS(1,1) == 'M'
Dataslash = '/';
elseif OS(1,1) == 'P'
Dataslash = '\';
end 
%% Processing Data 

% find matching pData files

%L3 Marvin's data
% f_Tm9LexA = find(i_FullField_OFF_3s_OFF_3s.*i_UASGCaMP6f_L3MH56Gal4_cross_to_w_mars.*~i_moving);

%L2 Marvins's data
f_Tm9LexA = find(i_FullField_OFF_3s_OFF_3s.*i_UASGCaMP6F_L221DhhGal4_cross_to_w_mars.*~i_moving);


% create structures by neuron
fTm9LexA = create_neuron_structure_all_ks(f_Tm9LexA);

% load all the data!
xTm9LexA = load_neuron_data10Hz_ks(fTm9LexA,pdatapath);

% This function groups the different interesting epoch combinations in a
% Matrix
mTm9LexA = aggregate_Adaptingcontrast(xTm9LexA,OFF_stimulus);
mTm9LexA.OFF_stimulus = OFF_stimulus;
%This function sorts the negative and positive correlated ROIs 
CorrelationA = crosscorrelation_negative(mTm9LexA); 
    
negativecorr = CorrelationA.negresponse;
    
negativecorr = logical(negativecorr);
Sample_size = length(unique(mTm9LexA.flyID));
    
DataSave = ('MH56_Control');
%% Test plots to check data
N_ROIs = size(mTm9LexA.rats,1);
N_EPOCHs = size(mTm9LexA.rats,2);
for cc = 1:N_EPOCHs
  rats(:,cc,:) = mTm9LexA.rats(negativecorr,cc,:);
  
end
close all;




cur_t5 = linspace(0,1,size(mTm9LexA.rats,2))';
cur_t = linspace(0,1,size(mTm9LexA.rats,3))';
AllTraces = figure('Units', 'pixels', 'Position',[100 100 500 375]); %Returns the graph with the subplots
color = varycolor(6); %This creates different color values for each epoch
for ii =1:N_EPOCHs-1

    hold on;
    epochRats = squeeze(rats(:,ii,:));
    [x,m,e] = mean_cat_full(epochRats,1,mTm9LexA.flyID); %m=mean; e=std
    h1 = plot_err_patch_v2((36*cur_t),m*100,e*100,color(ii,:),color(ii,:),'-',0); 
    %The last input parameter defines how transparent the error traces are, set 0
    %if no error traces are wanted
   
    if ii ~= N_EPOCHs
        set(gca, 'XTick', [])
    end
   hold on
end
%This makes the plot look nice
    ax = gca;
     ax.FontSize = 20;
      set(h1, 'linewidth'  ,2);
   set( gca,...
        'YGrid'       , 'on',...
        'XGrid'       , 'off',...
        'box', 'off',...
        'ylim', [-15 150],...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'off'      , ...
        'XColor'      , [.2 .2 .2], ...
        'YColor'      , [.2 .2 .2], ...
        'YTick'       , [0:10:250], ...
        'XTick'       , 0:5:35);
     xlabel('Time (s)');
    ylabel('dF/F - Calcium Signal');
    legend('hide');
    ylim ([0 20])
    ax.YLimMode = 'auto';
    %Save the plot to the destination set above under the name set above
% if autosave == 1
% print([num2str(Savelocation),Dataslash,'RawTraces_',num2str(DataSave),num2str(NDFilter)],'-depsc')
% end


%Here the epochs are plotted individually
for ii =1:N_EPOCHs
   figure(10*ii);
   hold on;
   
   epochRats = squeeze(rats(:,ii,120:230));
   [x,m,e] = mean_cat_full(epochRats,1,mTm9LexA.flyID); %m=mean; e=std
   h1 = plot_err_patch_v2((36*cur_t(120:230,1)),m,e,[0 0 1],[0.8 0.8 1]);
   
   xlabel('time (s)');
   ylabel('dF/F - Calcium Signal');
    ylim([-0.2 1]); 
   title(['Epoch ',num2str(ii)]);
   
   ax = gca;
   ax.FontSize = 20;
   legend('error',['n = ' num2str(Sample_size)]);
   grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , 0:0.2:2, ...
  'XTick'       , 10:5:25, ...
  'LineWidth'    , 1         );
  hold off;


end
%% Evaluating the decay rate to first step
%This section evaluates the decay rates of the A step of the stimulus as a mean over flies and ROIs 
if decay_analysis == 1;
fitresult=(1:N_EPOCHs-1);
gof=(1:N_EPOCHs-1);

for ii =1:N_EPOCHs-1
    cur_t = linspace(0,1,size(mTm9LexA.rats,3))';
epochRats = squeeze(rats(:,ii,:));
    [x,m,e] = mean_cat_full(epochRats,1,mTm9LexA.flyID); %m=mean; e=std
    %Define the time of the A step
    m = m(:,(150:180));
    e = e(:,(150:180));
    cur_t = cur_t((150:180),:);
    %Search for the max response and set it as the start of the traces
    [maxresponse,maxindex(ii)] = max (m);
    maxerror(ii) = e(maxindex(ii));
    m = m(:,(maxindex:end));
    e = e(:,(maxindex:end));
    cur_t = cur_t((maxindex:end),:);
    
    hold on
    %This saves the mean decay rate over ROIs for the 6 different steps
   [fitresult, gof,plotfit(:,ii)] = createFit(36*cur_t,m,DataSave,ndFilter);
   plotfit(1,ii).Color = color(ii,:);
   plotfit(1,ii).MarkerSize = 15;
   plotfit(2,ii).Color = 'k';
   decayrate(ii) = fitresult.b; %Here the actual data is saved
   
   
   %This saves the decay rates of different flies for the 6 different steps
     for ee = 1:length(x(:,1)) 
         clear xd
     cur_t = linspace(0,1,size(mTm9LexA.rats,3))';
    [x,m,e] = mean_cat_full(epochRats,1,mTm9LexA.flyID); %m=mean; e=std
    xd(1,:) = x(ee,(150:180));
    cur_t = cur_t((150:180),:);
    %Search for the maximal response and set it as the start of the trace
   [maxresponse_fly,maxindex_fly(ii,ee)] = max (xd(1,:));
    xd = xd(1,(maxindex_fly(ii,ee):end));
     cur_t = cur_t((maxindex_fly(ii,ee):end),:);
     
     if length(xd)> 10
    [fitresult, gof] = createFit(36*cur_t,xd(1,:),DataSave,ndFilter);
    decayrate_fly(ii,ee) = fitresult.b; %saves the decay rate for the individual flies
     end
   end
   %This plots the decay rate into the plot
text(16,0.5,['f(x)=a*exp(b*x) ','a = ',num2str(fitresult.a),' b = ', num2str(fitresult.b)]);
text(16,0,['adjusted R^2 = ', num2str(gof.adjrsquare)]);

ax = gca;
ax.FontSize = 20;
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , 0:0.2:2, ...
  'LineWidth'   , 1         );
%Save the plot
if autosave == 1
print([num2str(Savelocation),Dataslash,'Decay_RateGraph',num2str(DataSave),num2str(NDFilter),num2str(ii)],'-depsc')
end
hold off;

   
end
%Calculation the mean decay rate and the mean time of max response plus the
%std of these values
stdmaxindex = std(maxindex);
maxindex_avrg = mean(maxindex)/size(maxindex,2);
stddecayrate = std(decayrate);
meandecay = mean(decayrate);

close all
%Boxplot showing the time until max response for the 6 different epochs
Maxtimescat = figure('Units', 'pixels', 'Position',[100 100 500 375]);
hold on
boxplot(maxindex_fly');
 xlabel('Contrast (%)');
   ylabel('Time (s)');
    ylim([0 40]); 
   title('Time until max. response');
   
   ax = gca;
   ax.FontSize = 20;
   grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , 0:5:40, ...
  'XTick'       , 1:1:6, ...
  'LineWidth'    , 1         );
  hold off;
if autosave == 1
print([num2str(Savelocation),Dataslash,'Response_Time',num2str(DataSave),num2str(NDFilter),num2str(ii)],'-depsc')
end
%Boxplot showing the decar rates for the 6 different epochs
Mean_decayrate = figure('Units', 'pixels', 'Position',[100 100 500 375]);
hold on
boxplot(decayrate_fly');
 xlabel('Contrast (%)');
   ylabel('Decay rate (abritary units)');
    ylim([-2 2]); 
   title('Decay rates of responses');
   ax = gca;
   ax.FontSize = 20;
   grid on
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , -2:1:2, ...
  'XTick'       , 1:1:6, ...
  'LineWidth'    , 1         );
  hold off;
  
if autosave == 1
print([num2str(Savelocation),Dataslash,'Decay_Boxplot_Fly',num2str(DataSave),num2str(NDFilter),num2str(ii)],'-depsc')
end

end


%% Evaluating the responses to the first contrast step
 %here the response to the A steps of the stimulus will be evaluated. 
close all
clear ax
 THRESHOLD = 0;
fprintf('THRESHOLD for Signal sum: %d \n',THRESHOLD );

% take a mean of all assigned values, ignore 0
crossing_rats = rats; % Take only negative corelated ROIs
crossing_rats(crossing_rats==0) = NaN;

% Mean over all ROI's
mean_val = squeeze(nanmean(crossing_rats,1));

% Over fly ID
num_steps = size(mean_val,1);

mean_calciumStep = zeros(1,num_steps-1);
err_calciumStep = zeros(1,num_steps-1);


% Search max-Signal per ROI and take the mean of all max's per ROI signals
for kk=1:num_steps-1
    
    traces = squeeze(mTm9LexA.rats(negativecorr,kk,:)); %because mean_cat_full works only with 2 dimensions %taking all ROIs, the epoch we are one, full lenth (200)
    [x,m,e] = mean_cat_full(traces,1,mTm9LexA.flyID(negativecorr)); %m=mean; e=std
    
    %defining the time where we look for the baseline
    baseline = x(:,100:150); 
    
    %defining the time where we look for the response
    response = x(:,150:180); 

    decreasingCa = nanmean(baseline,2) - nanmean(response,2) >= 0; %nanmean is the mean of X, computed after removing NaN values
    increasingCa = nanmean(squeeze(baseline),2) - nanmean(squeeze(response),2) < 0;
    
    numFlies=size(x,1); %Added to prevent crash when number of flies does not equal number of Ca. traces. 
    calciumStep = zeros(numFlies,1);
    for ii = 1:numFlies
        if decreasingCa(ii)
            calciumStep(ii) = max(response(ii,:)) - mean(baseline(ii,:));  
        elseif increasingCa(ii) %this is necessary, if there are NaNs, both logicals are 0
            calciumStep(ii) = max(response(ii,:)) - mean(baseline(ii,:));  %        
        else
            calciumStep(ii) = NaN;
        end
    end
    
    % take the mean off all steps
    mean_calciumStep(kk) = mean(calciumStep);
    %This saves the data for individual flies for later ANOVa analysis
    Calcium_stepforAnova(:,kk) = calciumStep;
    Calcium_stepforAnova_test{kk} = calciumStep;
    % Standard error
    err_calciumStep(kk) = nanstd(calciumStep)/sqrt(length(calciumStep));
    %Save the Data for analysis
    
end

%% Contrast step Over fly


num_steps = size(mean_val,1);

numFlies = length(unique(mTm9LexA.flyID));

mean_ContrastStep = zeros(1,num_steps-1);
err_ContrastStep = zeros(1,num_steps-1);



% Search max-Signal per ROI and take the mean of all max's per ROI signals
for kk=1:num_steps-1
    
    traces = squeeze(mTm9LexA.rats(negativecorr,kk,:)); %because mean_cat_full works only with 2 dimensions %taking all ROIs, the epoch we are one, full lenth (200)
    [x,m,e] = mean_cat_full(traces,1,mTm9LexA.flyID(negativecorr)); %m=mean; e=std
    
    %ks: defining the time where we look for the baseline
    baseline = x(:,100:150); 

    %ks: defining the time where we look for the response
    response = x(:,180:210); 
  
    decreasingCa = nanmean(baseline,2) - nanmean(response,2) >= 0; %nanmean is the mean of X, computed after removing NaN values
    increasingCa = nanmean(squeeze(baseline),2) - nanmean(squeeze(response),2) < 0;
    
    numFlies=size(x,1);%mars Added to prevent crash when number of flies does not equal number of Ca. traces. 
    ContrastStep = zeros(numFlies,1);
    for ii = 1:numFlies
       %this is necessary, if there are NaNs, both logicals are 0
            ContrastStep(ii) = max(response(ii,:)) - mean(baseline(ii,:));  %        
       
    end
    
    all_step_samples{kk} = ContrastStep; %ks for ANOVA
    
    % take the mean off all steps
   mean_ContrastStep(kk) = nanmean(ContrastStep);
    err_ContrastStep(kk) = nanstd(ContrastStep)/sqrt(length(calciumStep));
     
end

%% Last ON  step in the OFF_OFF stimulus
 %Here the amplitude of the response to the ON step after the B step is
 %calculated

close all
 THRESHOLD = 0;
fprintf('THRESHOLD for Signal sum: %d \n',THRESHOLD );

% take a mean of all assigned values, ignore 0
crossing_rats = rats; % Take only negative corelated ROIs
crossing_rats(crossing_rats==0) = NaN;

% Mean over all ROI's
mean_val = squeeze(nanmean(crossing_rats,1));

numROIs = size(crossing_rats,1);
num_steps = (size(mean_val,1));


mean_LastStep = zeros(1,num_steps);
err_LastStep = zeros(1,num_steps);

% Search max-Signal per ROI and take the mean of all max's per ROI signals
for kk=1:num_steps
    
    %defining the time where we look for the baseline
    baseline_Contrast = crossing_rats(:,kk,(120:150));
    baseline_Contrast = squeeze(baseline_Contrast);
    %defining the time where we look for the response
    response_Contrast = crossing_rats(:,kk,(210:230));
    response_Contrast = squeeze(response_Contrast);
    
    
    LastStep = zeros(numROIs,1);
    for ii = 1:numROIs
        %this is necessary, if there are NaNs, both logicals are 0
            LastStep(ii) = min(response_Contrast(ii,:)) - mean(baseline_Contrast(ii,:));  %        
       
    end
    
    % take the mean of all steps
    mean_LastStep(kk) = nanmean(LastStep);
    % Standard error
    err_LastStep(kk) = nanstd(LastStep)/sqrt(length(LastStep));
     
    
end


%% B step calculated from baseline A step

close all
 THRESHOLD = 0;
fprintf('THRESHOLD for Signal sum: %d \n',THRESHOLD );

% take a mean of all assigned values, ignore 0
%crossing_rats = mTm9LexA.rats; % Take all ROIs
crossing_rats = rats; % Take only negative corelated ROIs
crossing_rats(crossing_rats==0) = NaN;

% Mean over all ROI's
mean_val = squeeze(nanmean(crossing_rats,1));

numROIs = size(crossing_rats,1);
num_steps = (size(mean_val,1));

mean_BStep = zeros(1,num_steps);
err_BStep = zeros(1,num_steps);

% Search max-Signal per ROI and take the mean of all max's per ROI signals
for kk=1:num_steps
    
    %ks: defining the time where we look for the baseline
    baseline_BStep = crossing_rats(:,kk,(170:190));
    baseline_BStep = squeeze(baseline_BStep);
    %ks: defining the time where we look for the response
    response_BStep = crossing_rats(:,kk,(180:210));
    response_BStep = squeeze(response_BStep);
    
    
    BStep = zeros(numROIs,1);
    for ii = 1:numROIs
        %this is necessary, if there are NaNs, both logicals are 0
               BStep(ii) = max(response_BStep(ii,:)) - min(baseline_BStep(ii,:));  %        
       
    end
    
    % take the mean of all steps
    mean_BStep(kk) = nanmean(BStep);
    % Standard error
    err_BStep(kk) = nanstd(BStep)/sqrt(length(BStep));
     
    
end
%% Graphs

%Here all the data is plotted in several different ways. 
y1 = mean_calciumStep; %ks
x1 = real_luminance; %Difference_luminance;  %real_luminance;
e1 = err_calciumStep;


%Plot the A Step response over luminance 
%The informations for the legend will be returned in Luminanceplot
figure('Units', 'pixels', 'Position',[100 100 500 375]);
grid on;
hold on;
Luminanceplot(1) = scatter(x1,y1,75,'filled','MarkerFaceColor',[0.9 0.3 0.2],'DisplayName','Response to A step');
Luminanceplot(2) = errorbar(x1,y1,e1,'color',[0.9 0.3 0.2],'Linestyle','none','DisplayName','SEM');

%Plot a regression line
[coeffs1,S1] = polyfit(x1, y1, 1);
%Calculate R square
f1L = polyval(coeffs1,x1);
[r1L, rmse1L] = rsquare(y1,f1L);
% Get fitted values
fittedX1 = linspace(min(x1), max(x1), 100);
fittedY1 = polyval(coeffs1, fittedX1);
% Plot the fitted line
hold on;
Luminanceplot(3) =  plot(fittedX1, fittedY1, 'k-.', 'LineWidth', 2,'DisplayName','Linear regression');
Error_Fit1 = polyval(coeffs1,x1);
T1 = table(x1,y1,Error_Fit1,y1-Error_Fit1,'VariableNames',{'X','Y','Fit','FitError'});
% Calculate the confidence values and plot the confidence lines 
% The error is not yet included
[Y1,DELTA1] = polyconf(coeffs1,fittedX1,S1,'alpha',0.05);
% fittedX1 = sort(fittedX1,'descend');
Luminanceplot(4) = plot(fittedX1,Y1+DELTA1,'--','color',[0.9 0.3 0.2],'DisplayName','95% confidence interval');
plot(fittedX1,Y1-DELTA1,'--','color',[0.9 0.3 0.2]);


%Plot the response of the B step
%Change of luminance depending on the ND filters
y2 = mean_ContrastStep; %ks
x2 =  real_luminance1; % Difference_luminance1;  %real_luminance1;
e2 = err_ContrastStep;


hold on
Luminanceplot(5) = scatter(x2,y2,80,'s', 'filled','MarkerFaceColor',[0.1 0.6 0.8] ,'DisplayName','Response to B step'); %square
hold on
Luminanceplot(6) =  errorbar(x2,y2,e2,'color',[0.1 0.6 0.8],'Marker','none','Linestyle', 'none','MarkerFaceColor','b','DisplayName','SEM');
hold on

[coeffs2,S2] = polyfit(x2, y2, 1);
% Get fitted values
fittedX2 = linspace(min(x2), max(x2), 100);
fittedY2 = polyval(coeffs2, fittedX2);
%Calculate R sqaure
f2L = polyval(coeffs2,x2);
[r2L, rmse2L] = rsquare(y2,f2L); 
% Plot the fitted line
hold on;
Luminanceplot(7) = plot(fittedX2, fittedY2, 'k--', 'LineWidth', 2,'DisplayName','Linear regression');
%Calculate the confidence values and plot the confidence lines 
%The error is not yet included
[Y2,DELTA2] = polyconf(coeffs2,fittedX2,S2,'alpha',0.05);
Luminanceplot(8) = plot(fittedX2,Y2+DELTA2,'--','color',[0.1 0.6 0.8],'DisplayName','95 % confidence interval');
plot(fittedX2,Y2-DELTA2,'--','color',[0.1 0.6 0.8], 'DisplayName','none');



hold on 
xlabel('Luminance (%)');
fprintf(['n = ' num2str(Sample_size)]);
title([num2str(DataSave),' response to A and B steps over luminance']);
ylabel('dF/F');
legend(Luminanceplot(1:8),'location','best');
r1L = num2str(r1L);
r2L = num2str(r2L);
%This prints the equation values on the plot
if r1L ~= 0
    text(fittedX1(100)-7,fittedY1(100)+0.15,['R^{2}= ',r1L(1:4)],'fontsize',12);
    text(fittedX1(100)-7,fittedY1(100)+0.1,['y =', num2str(coeffs1(1)),'*x+',num2str(coeffs1(2))],'fontsize',12);
end
if r2L ~=0 
    text(fittedX2(100)-2,fittedY2(100)-0.2,['R^{2}= ',r2L(1:4)],'fontsize',12);
    text(fittedX2(100)-2,fittedY2(100)-0.25,['y =', num2str(coeffs2(1)),'*x+',num2str(coeffs2(2))],'fontsize',12);
end

ylim([-0.2 1.0]);
% ylim([-0.6 2.5]);
% ylim([-0.4 3.0]);
ax = gca;
ax.FontSize = 20;
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , 0:0.2:1.0, ...
  'LineWidth'   , 1         );
if autosave == 1
print([num2str(Savelocation),Dataslash,'Response_over_Luminance_',num2str(DataSave),num2str(NDFilter)],'-depsc')
end
hold off;

%% Plot ALL steps over contrast

figure('Units', 'pixels', 'Position',[100 100 500 375]);
grid on;

% Plot the data only for OFF steps
y1 = mean_calciumStep; %ks
x1 = real_contrast1;
e1 = err_calciumStep;

%The data for the legend is returned in ContrastGra
hold on;
ContrastGra(1) = scatter(x1,y1,80,'filled','MarkerFaceColor',[0.9 0.3 0.2],'DisplayName','Response to A step');
ContrastGra (2) = errorbar(x1,y1,e1,'color',[0.9 0.3 0.2],'Marker','none','Linestyle','none','DisplayName','SEM');

%Plot a regression line
[coeffs1,S1] = polyfit(x1, y1, 1);
%Calculate R
f1c = polyval(coeffs1,x1);
[r1c, rmse1c] = rsquare(y1,f1c);

% Get fitted values
fittedX1 = linspace(min(x1), max(x1), 100);
fittedY1 = polyval(coeffs1, fittedX1);
% Plot the fitted line
hold on;
ContrastGra(3) = plot(fittedX1, fittedY1, 'k-.', 'LineWidth', 2,'DisplayName','Linear regression');
Error_Fit1 = polyval(coeffs1,x1);
T1 = table(x1,y1,Error_Fit1,y1-Error_Fit1,'VariableNames',{'X','Y','Fit','FitError'});
% Calculate the confidence values and plot the confidence lines 
% The error is not yet included
[Y1,DELTA1] = polyconf(coeffs1,fittedX1,S1,'alpha',0.05);
ContrastGra(4) = plot(fittedX1,Y1+DELTA1,'--','color',[0.9 0.3 0.2],'DisplayName', '95% confidence interval');
plot(fittedX1,Y1-DELTA1,'--','color',[0.9 0.3 0.2]);



% Plot the data only for ON steps
y2 = mean_ContrastStep; %ks
x2 = real_contrast;
e2 = err_ContrastStep;

ContrastGra (5) = scatter(x2,y2,80,'s','filled','MarkerFaceColor',[0.1 0.6 0.8],'DisplayName','Response to B step');
ContrastGra (6) = errorbar(x2,y2,e2,'color',[0.1 0.6 0.8],'Marker','none','Linestyle','none','DisplayName','SEM');

[coeffs2,S2] = polyfit(x2, y2, 1);
% Get fitted values
fittedX2 = linspace(min(x2), max(x2), 100);
fittedY2 = polyval(coeffs2, fittedX2);
f2c = polyval(coeffs2,x1);
[r2c, rmse2c] = rsquare(y2,f2c); 
% Plot the fitted line
hold on;
ContrastGra(7) = plot(fittedX2, fittedY2, 'k--', 'LineWidth', 2,'Displayname', 'Linear regression');
%Calculate the confidence values and plot the confidence lines 
%The error is not yet included
[Y2,DELTA2] = polyconf(coeffs2,fittedX2,S2,'alpha',0.05);
ContrastGra(8) = plot(fittedX2,Y2+DELTA2,'--','color',[0.1 0.6 0.8], 'DisplayName','95% confidence interval');
plot(fittedX2,Y2-DELTA2,'--','color',[0.1 0.6 0.8]);

Error_Fit2 = polyval(coeffs2,x2);
T2 = table(x2,y2,Error_Fit2,y2-Error_Fit2,'VariableNames',{'X','Y','Fit','FitError'});
hold on 
% ylim([-0.9 1.5]);
hold off
title([num2str(DataSave),' response to linear increasing and 25% contrast steps']);

% Add titles

xlabel('Contrast (%)');
ylim([-0.2 1.0]);
% ylim([-0.6 2.5]);
% ylim([-0.4 3.0]);
xlim([-90 -20]);
ylabel('dF/F - Calcium Signal');
legend (ContrastGra(1:8),'location','best');
r1L = num2str(r1L);
r2L = num2str(r2L);
if r1c ~= 0
    text(fittedX1(1)+2,fittedY1(1)+0.1,['R^{2}= ',r1L(1:4)],'fontsize',12);
    text(fittedX1(1)+2,fittedY1(1)+0.05,['y =', num2str(coeffs1(1)),'*x+',num2str(coeffs1(2))],'fontsize',12);
end
if r2c ~=0 
    text(fittedX2(1)+2,fittedY2(1)-0.1,['R^{2}= ',r2L(1:4)],'fontsize',12);
    text(fittedX2(1)+2,fittedY2(1)-0.15,['y =', num2str(coeffs2(1)),'*x+',num2str(coeffs2(2))],'fontsize',12);
end
ax=gca;
ax.FontSize = 20;
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , 0:0.2:1.0, ...
  'LineWidth'   , 1         );
if autosave ==1
print([num2str(Savelocation),Dataslash,'Response_over_Contrast_',num2str(DataSave),num2str(NDFilter)],'-depsc')
end
hold off;



%% ANOVA test - A step
%A step data
Calcium_stepforAnova = Calcium_stepforAnova_test;
anovaVector = [];
groupVector = [];

for iAnovaPrep = 1:length(Calcium_stepforAnova)
    anovaVector = [anovaVector Calcium_stepforAnova{iAnovaPrep}'];
    groupVector = [groupVector zeros(length(Calcium_stepforAnova{iAnovaPrep}),1)'+iAnovaPrep];    
end

%prepare for ANOVA

[p,tbl,stats] = anova1(anovaVector,groupVector);

pause
[results,means] = multcompare(stats,'CType','bonferroni')
% [results,means] = multcompare(stats,'CType','dunn-sidak')
results(end)



%% ANOVA test - B step
%B step data
Contrast_step_Bstep = all_step_samples;
anovaVector = [];
groupVector = [];

for iAnovaPrep = 1:length(Contrast_step_Bstep)
    anovaVector = [anovaVector Contrast_step_Bstep{iAnovaPrep}'];
    groupVector = [groupVector zeros(length(Contrast_step_Bstep{iAnovaPrep}),1)'+iAnovaPrep];
end

%prepare for ANOVA

[p,tbl,stats] = anova1(anovaVector,groupVector);

pause
[results,means] = multcompare(stats,'CType','bonferroni')
% [results,means] = multcompare(stats,'CType','dunn-sidak')
results(end)



