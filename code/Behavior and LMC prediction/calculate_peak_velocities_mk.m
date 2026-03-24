function []=calculate_peak_velocities_mk(flydata,outputfile)
%flydata is a 1Xn struct array with n=#flies.
%length(flydata(1).acc_dtheta) is the #groups (number of different OFF edge categories)
%groupMeanperFly are peak velocities of the trial-averaged time trace per fly per OFF edge category.
%groupMean and groupSEM are statistics of the 'groupMeanperFly' data,
%after baseline correction, calculated over flies.
%baseline is calculated per fly per group from the last 0.2 s of preceeding inter-stimulus intervals. 

%% constants
motionStart=0.7; %in s. This is because each row of the data, corresponding to one 
                 % stimulus instance,includes 0.2 s before the instance started, 
                 % followed by 0.5 s static background presentation. 
                 % Thus, OFF edge motion starts at 0.2+0.5=0.7 s.
durMotion=0.75; %in s.
durBaseline=0.2; % the 0.2 s preceeding each stimulus instance for baseline correction.

%% processing
groupMeanperFly=zeros(size(flydata,2),length(flydata(1).acc_dtheta));
baselineVels=zeros(size(flydata,2),length(flydata(1).acc_dtheta));
for i=1:length(flydata(1).acc_dtheta) % number of groups
    for j=1:size(flydata,2) %number of flies
        baselineVels(j,i)=max(mean(flydata(j).acc_dtheta{i}...
            (:,1:round(durBaseline*120)),1));
        groupMeanperFly(j,i)=max(mean(flydata(j).acc_dtheta{i}...
            (:,round((0.1+motionStart)*120):round((0.1+motionStart+durMotion)*120)),1));
    end
end

groupMeanperFlycorrected=groupMeanperFly-baselineVels; %baseline correction
groupMean=mean(groupMeanperFlycorrected); %mean across all flies
groupSEM=std(groupMeanperFlycorrected)/sqrt(length(flydata));

%% saving
save(['/your_peak_velocity_datapath/',outputfile,'.mat'],...
    'groupMean','groupSEM','groupMeanperFly','groupMeanperFlyCorrected',...
    'peakVelsperFlyperGroup','baselineVels');

end

    
