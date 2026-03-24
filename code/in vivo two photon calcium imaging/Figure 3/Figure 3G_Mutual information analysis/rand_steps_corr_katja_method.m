function negative_correlated = rand_steps_corr_katja_method(neuronRawData)
% From Katja Sporar's negative correlation method used in Ketkar_Sporar
% for this stimulus.
aggregated_data = aggregate_05StepsON_OFF(neuronRawData); %taking grey as a baseline

numROIs = size(aggregated_data.rats,1);


% There are 20 epochs which go from 0->1, 0->2 ...
from = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4];
to = [1,2,3,4,0,2,3,4,0,1,3,4,0,1,2,4,0,1,2,3];

step = to-from;
step(step<0) = 0; 
on_step = logical(step); 

negative_correlated = zeros(numROIs,1);
positive_correlated = zeros(numROIs,1);

% Check which cells are positively or negatively correlated
for ii=1:numROIs
    responses_roi = squeeze(aggregated_data.rats(ii,:,:));
    
    % Mean signal over all On Steps or Off Steps
    mean_onStep = nanmean(responses_roi(on_step,:),1);
    mean_offStep = nanmean(responses_roi(~on_step,:),1);
   
    reaction_onStep = nanmean(mean_onStep(100:120))-nanmean(mean_onStep(80:100)); % if <0 then neg. correlated, if >0 then pos. correlated
    reaction_offStep = nanmean(mean_offStep(100:120))-nanmean(mean_offStep(80:100)); % if >0 then neg. correlated, if <0 then pos. correlated 
    
    if reaction_offStep < 0 && reaction_onStep > 0
       % positive correlated
       positive_correlated(ii)=1;
    elseif reaction_offStep > 0 && reaction_onStep < 0
        % negative correlated
        negative_correlated(ii)=1;
    else 
        % weird stuff, better drop this ROI
    end

end

% Array stores 'true' (1) if this ROI is negatively correlated
negative_correlated = logical(negative_correlated);
