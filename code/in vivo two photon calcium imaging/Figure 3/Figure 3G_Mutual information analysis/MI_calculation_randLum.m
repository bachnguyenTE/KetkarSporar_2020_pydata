

function neuronDataOfInterest = Sporar_et_al_MI_calculation_randLum(neuronRawData)
% This function calculates mutual information and correlation 
% Used for Sporar et al.

%% For each ROI MI separately
FPS = 10;
iGoodROI = 1;
all_flyIDs = [neuronRawData.flyID];


negCorrROIs = neuronRawData(1).negCorrROIs';

% Nearest neighbor estimation parameter K
K = 20; % manual entry
for iROI=1:length(neuronRawData)
    
    % Getting the raw trace, not using the first 5s
    rawTrace = neuronRawData(iROI).iratio;
    rawTrace = rawTrace(5*FPS:end); % getting rid of first 5s (possible noise)
 
    % dF/F calculation with the mean of the trace
    relativeTrace = rawTrace./mean(rawTrace) - 1; 
     
    % Smoothing with a moving average of 5 to get rid of some fast noise
    relativeTrace = smooth(relativeTrace);
    relativeTrace = relativeTrace';
    
    % Absolute thresholding and negative correlation thresholding
    if  abs(max(relativeTrace(10*FPS:end)))<0.3 
        continue
    elseif negCorrROIs(iROI)==0
        continue     
    end
    
    % Stimuli processing
    stimulus = neuronRawData(iROI).istim(5*FPS:end);
    luminanceTrace = stimulus; 

    % Scaling responses for MI Toolbox
    scaling_factor_lum = length(unique(luminanceTrace));
    scaled_trace = round(relativeTrace * scaling_factor_lum);
    
    %%%%% Random stimulus MI calculation %%%%% 
    
    % Two randomization methods:
    % 1) Only randomizing the 10s epochs together keeping the temporal
    % structure in the stimulus.
    % 2) Randomizing all values and getting rid of the temporal structure in
    % the stimulus
    for iRand = 1:100 %100 times randomizing
        
        % Keeping temporal structure
        xLum1=randperm(numel(luminanceTrace));
        shuffLumVals1 =luminanceTrace(xLum1);
        shuffLumVals1 = [0 shuffLumVals1 0];
        changeIndx = find(diff(luminanceTrace));
        startIdx = 1;
        shuffLumTrace1 = luminanceTrace;
        for iShuffle = 1:numel(changeIndx)+1
            if iShuffle > numel(changeIndx)
                shuffLumTrace1(changeIndx(iShuffle-1):end) = shuffLumVals1(iShuffle);
            else
                shuffLumTrace1(startIdx:changeIndx(iShuffle)) = shuffLumVals1(iShuffle);
                startIdx = changeIndx(iShuffle);
            end
        end
        % Stimulus values randomization
        Xlum2=randperm(numel(luminanceTrace));
        
        % All randomization
        shuffledLuminanceTrace2 =luminanceTrace(Xlum2);
        
        randLumMatrix_dis_cont1(iRand) = discrete_continuous_info_fast(shuffLumTrace1,relativeTrace,K);
        randLumMatrix_MI_toolbox1(iRand) = mi(scaled_trace',shuffLumTrace1');
%         
%         randLumMatrix_dis_cont2(iRand) = discrete_continuous_info_fast(shuffledLuminanceTrace2,relativeTrace,K);
%         randLumMatrix_MI_toolbox2(iRand) = mi(scaled_trace',shuffledLuminanceTrace2');
%         
%         randLumMatrix_dis_cont1(iRand) = 0;
%         randLumMatrix_MI_toolbox1(iRand) = 0;
%       % First method used so just assign 0 to second
        randLumMatrix_dis_cont2(iRand) = 0;
        randLumMatrix_MI_toolbox2(iRand) = 0;

    end

    miLumTraceRandomized_MI_toolbox1(iGoodROI) = mean(randLumMatrix_MI_toolbox1);
    miLumTraceRandomized_dis_cont1(iGoodROI) = mean(randLumMatrix_dis_cont1);

    miLumTraceRandomized_MI_toolbox2(iGoodROI) = mean(randLumMatrix_MI_toolbox2);
    miLumTraceRandomized_dis_cont2(iGoodROI) = mean(randLumMatrix_dis_cont2);

    %%%%% Mutual information calculations %%%%% 
    
    % MI Toolbox Brown et al. 2012, scaled responses are used
    miLuminance_whole_trace_scaled_MI1(iGoodROI) = mi( scaled_trace', luminanceTrace' );
   
    % MI Brian C. Ross 2014, no scaling needed
    miLuminance_whole_trace_MI2(iGoodROI) = discrete_continuous_info_fast...
        (luminanceTrace,relativeTrace,K);
   
    % Saving the rest of the information from the data
    corrLuminance(iGoodROI) = corr( scaled_trace', luminanceTrace' );
    flyIDs(iGoodROI) = all_flyIDs(iROI);
    iGoodROI = iGoodROI + 1;
   
end   
fprintf('Pass rate of current genotype: %.2f', (iGoodROI-1)/iROI)

%%


neuronDataOfInterest.flyID = flyIDs;
neuronDataOfInterest.corr = corrLuminance ;

neuronDataOfInterest.Ilum1 = miLuminance_whole_trace_scaled_MI1;
neuronDataOfInterest.Ilum2 = miLuminance_whole_trace_MI2;

% neuronDataOfInterest.Ilum_early1 = miLuminance_early_trace_lum_MI1;
% neuronDataOfInterest.Ilum_early2 = miLuminance_early_trace_MI2;

neuronDataOfInterest.Ilum1_rand1 = miLumTraceRandomized_MI_toolbox1;
neuronDataOfInterest.Ilum2_rand1 = miLumTraceRandomized_dis_cont1;

neuronDataOfInterest.Ilum1_rand2 = miLumTraceRandomized_MI_toolbox2;
neuronDataOfInterest.Ilum2_rand2 = miLumTraceRandomized_dis_cont2;

