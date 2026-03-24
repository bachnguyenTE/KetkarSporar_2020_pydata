

function neuronDataOfInterest = MI_time_analysis(neuronRawData)
% This function calculates mutual information 

%% For each ROI MI separately
FPS = 10;
iGoodROI = 1;
all_flyIDs = [neuronRawData.flyID];


negCorrROIs = neuronRawData(1).negCorrROIs';

% Nearest neighbor estimation parameter K
K = 20; % manual entry
 for iROI=1:length(neuronRawData)
   
    % Getting the raw trace, not using the first 50 frames
    rawTrace = neuronRawData(iROI).iratio;
    rawTrace = rawTrace(5*FPS:end); % getting rid of first 5s
 
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
    ds = diff(luminanceTrace); % finding stimulus changes
    b_inds = find(ds)+1; % beginning of each distinct stimuli
    
    % Doing time analysis with window of WINDOWSIZE * 100ms
    window_size = 2*FPS; % in frames 2 second windows %
    step_size = 5; % Going in 0.5s steps 
    curr_window = 0;
    startIndices = 1:step_size:81;
    
   
    
    for iWindow = 1:length(startIndices)
        % Creating traces 
        timeAnalysisRespTrace = zeros(length(b_inds)-1, window_size);
        timeAnalysisLumTrace = zeros(length(b_inds)-1, window_size);
        timeAnalysisContTrace = zeros(length(b_inds)-1, window_size);
        startIndx = startIndices(iWindow);
        endIndx = startIndx +  window_size-1;
        
        for iEpochs = 1:(length(b_inds)-1)
            
            % Taking the current epochs luminance and response traces (put -5
            % since I want to be sure that its always within the epoch
            currT = relativeTrace(b_inds(iEpochs):...
             b_inds(iEpochs+1));
            currLuminanceTrace = luminanceTrace(b_inds(iEpochs):...
             b_inds(iEpochs+1));
         
            timeAnalysisRespTrace(iEpochs,:) = currT(startIndx : endIndx);
            timeAnalysisLumTrace(iEpochs,:) = currLuminanceTrace(startIndx : endIndx);
        end
        curr_window = curr_window + 1;
        totalLength = size(timeAnalysisRespTrace,1) * size(timeAnalysisRespTrace,2);
        timeAnalysisRespTrace = reshape(timeAnalysisRespTrace',totalLength,1);
        timeAnalysisLumTrace = reshape(timeAnalysisLumTrace',totalLength,1);
        
        
        % Scaling for MI Toolbox
        scaling_factor_lum = length(unique(luminanceTrace));
        timeAnalysisRespTrace_luminance_scaled = round((timeAnalysisRespTrace./max(relativeTrace)) * scaling_factor_lum);
        
        % Randomizing LUMINANCE
        % all numbers
        Xlum=randperm(numel(timeAnalysisLumTrace));
        shuffledLuminanceTrace1 =reshape(timeAnalysisLumTrace(Xlum),size(timeAnalysisLumTrace));    
        
        % Only epoch randomization
        xLum=randperm(numel(timeAnalysisLumTrace));
        shuffLumVals1 =reshape(timeAnalysisLumTrace(xLum),size(timeAnalysisLumTrace));
        shuffLumVals1 = [0 shuffLumVals1' 0];
        changeIndx = find(diff(timeAnalysisLumTrace));
        startIdx = 1;
        shuffLumTrace = timeAnalysisLumTrace;
        for iShuffle = 1:numel(changeIndx)+1
            if iShuffle > numel(changeIndx)
                shuffLumTrace(changeIndx(iShuffle-1):end) = shuffLumVals1(iShuffle);
            else

                shuffLumTrace(startIdx:changeIndx(iShuffle)) = shuffLumVals1(iShuffle);
                startIdx = changeIndx(iShuffle);
            end


        end
        
   
    
        
        
        
        % MI Toolbox method
        MItimeTraceLumMI1(iGoodROI, iWindow) = mi...
            (timeAnalysisRespTrace_luminance_scaled, timeAnalysisLumTrace);
        MItimeTraceLumRandMI1(iGoodROI, iWindow) = mi...
            (timeAnalysisRespTrace_luminance_scaled, shuffLumTrace);
        MItimeTraceLumRand2MI1(iGoodROI, iWindow) = mi...
            (timeAnalysisRespTrace_luminance_scaled, shuffledLuminanceTrace1);
        
        % MI Brian C. Ross 2014, no scaling needed
         
        % Luminance
        MItimeTraceLumMI2(iGoodROI, iWindow) = discrete_continuous_info_fast...
            (timeAnalysisLumTrace,timeAnalysisRespTrace,K);
        MItimeTraceLumRandMI2(iGoodROI, iWindow) = discrete_continuous_info_fast...
            ( shuffLumTrace,timeAnalysisRespTrace,K );
        MItimeTraceLumRand2MI2(iGoodROI, iWindow) = discrete_continuous_info_fast...
            ( shuffledLuminanceTrace1,timeAnalysisRespTrace,K );

    end
    flyIDs(iGoodROI) = all_flyIDs(iROI);
    iGoodROI = iGoodROI + 1;
end   

neuronDataOfInterest.flyID = flyIDs;
neuronDataOfInterest.MItimeTraceLum = MItimeTraceLumMI2 ;
neuronDataOfInterest.MItimeTraceLumRand = MItimeTraceLumRandMI2 ;

neuronDataOfInterest.MItimeTraceLumMI_Toolbox = MItimeTraceLumMI1 ;
neuronDataOfInterest.MItimeTraceLumMI_Toolbox_rand = MItimeTraceLumRandMI1 ;





