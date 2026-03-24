function alignSingleTimeSeries_ms(dirPath, wantMovie)
    
    if nargin < 2
        wantMovie = false;
    end
    
    cd(dirPath)
    out.fileloc = dirPath;
    [~, flyID] = fileparts(fileparts(dirPath));
    out.dataID = flyID;

    %%

    find the .lsm file  %%%%%%%% We are not using any .lsm file in the paper %%%%%%%%%%%
    d=dir('*.lsm');
    filename = [d.name];

    % If there is no .lsm file look for .xml file.
    if isempty(filename)
        folder = dir('*.xml');
        filename = [folder.name];
        if isempty(filename)
            xmlDir = dir('T*');
            xmlDirName = xmlDir([xmlDir.isdir]).name;
            cd(xmlDirName)
            folder = dir('*.xml');
            filename = [folder.name];
            filetype = 'xml';
        else
            filetype = 'xml';
        end
    else    
        filetype = 'lsm';
    end
    
    %% Extract information from imaging file

    switch filetype
        case 'lsm'
            [lsminf, scanInfo, ~] = lsminfo(filename);
            imagingInfo = lsminf;
        case 'xml'
            [imagingInfo, scanInfo] = getXmlInfo(filename); %information for every single frame, as well as scan settings
    end

    %% Information of stimulus time series.

    % Number of frames.
    nImages = numel(imagingInfo.stimulusFrames); 
    % Image depth in bits.
    imageDepth = str2double(scanInfo.bitDepth);
    % Image size in number of pixels.
    height = str2double(scanInfo.linesPerFrame);
    width = str2double(scanInfo.pixelsPerLine);

    %% Define output file.
    % Remove extension from filename.
    nameNoExt = regexp(filename, ['.' filetype], 'split');
    % Create filename for aligned images, with TIFF format.
    alignedFile = [nameNoExt{1}, '_aligned.tif'];
    % Save time stamps.
    timeStampsFile = [nameNoExt{1} '_times.mat'];
    % Extract the time stamps.
    stimulusFrames = [imagingInfo.stimulusFrames{:}];
    % Compute mean frame duration from relative frame times. 
    relativeFrameTimes = {stimulusFrames(:).relativeTime};
    relativeFrameTimes = str2double(relativeFrameTimes');
    relFrameLength = mean(diff(relativeFrameTimes)); %%%%%%%%%%%%ks: we don' use this - DELETE???%%%%%%%%%%%
    % Compute mean frame duration from relative frame times. 
    absoluteFrameTimes = {stimulusFrames(:).absoluteTime};
    absoluteFrameTimes = str2double(absoluteFrameTimes');

    absFrameLength = mean(diff(absoluteFrameTimes)); %%%%%%%%%%%%ks: we don' use this - DELETE???%%%%%%%%%%%
    % Frame duration from microscope settings.
    framePeriod = str2double(scanInfo.framePeriod);
    % Time stamps, relative or absolute?.
    ts = relativeFrameTimes;
    % ts = imagingInfo.TimeStamps.TimeStamps;
    save(timeStampsFile, 'ts');
    % fr = mean(diff(ts));

    %% Read tiff time series.

    % Read the reference frames.
    switch filetype
        case 'lsm'
            imageArray = tiffread(filename,1: 2: 2 * nImages);
        case 'xml'
            imageArray = readTwoPhotonTimeSeries(filename, imagingInfo);
        otherwise
            error(['Unrecognized file type: ' filetype])
    end

    %% Perform alignment.
    % Create reference stack, using the first 30 frames.
    refStack = zeros(height, width, 30);
    switch filetype
        case 'lsm'
            for iFrame = 1: 30
                refStack(:, :, iFrame) = imageArray(iFrame).data;
            end
            % Maximum intensity projection of reference stack.
            refFrame = max(refStack, [], 3);
            ref = uint8(refFrame);
            % Read reference image.
            im = tiffread(filename, refFrame*2-1);
            refFrame = im2double(im.data);
        case 'xml'
            for iFrame = 1: 30
                refStack(:, :, iFrame) = imageArray(:, :, iFrame);
            end
            % Maximum intensity projection of reference stack.
            refFrame = max(refStack, [], 3);
        otherwise
            error(['Unrecognized file type: ' filetype])        
    end
    %% Perform alignment based on maximizing image crosscorrelation in the Fourier space.
    in1 = imageArray;
    in2 = refFrame;
    in3 = filetype;
    in4 = alignedFile;
    
    [out1, out2, out3] = fourierCrossCorrelAlignment(in1, in2, in3, in4);
    
    registeredImages = out1;
    unregisteredImages = out2;
    medianFilteredImages = out3;
    
    % Store unaligned image in output.
    out.ch1 = im2double(unregisteredImages);
    % Store aligned image in output.
    out.ch1a = im2double(registeredImages);
    % Store aligned filtered image in output.
    out.ch1b = medianFilteredImages;

    %% Write to out.xml the info in scanInfo struct.
    out.xml = scanInfo;
    out.xml.frames = nImages;
    % Assign frametime, frame rate, and z depth for subsequent scripts.
    out.xml.frametime = framePeriod;
    out.xml.framerate = 1/framePeriod; 
%     out.xml.zdepth = str2double(scanInfo.positionCurrent.ZAxis); 
    out.xml.zdepth = str2double(scanInfo.positionCurrent.Z); %after z piezo installation   
    out.xml.absoluteFrameTimes = absoluteFrameTimes;

    %% Process stimulus file.
    cd ..
    d = dir('_stimulus_output*');
    A = importdata(d.name);
    ch3 = A.data(:,4); 
    frame_nums = A.data(:,8); 

    fid=fopen(d.name,'r');
    currline=fgetl(fid);
    ind = strfind(currline,'\');
    stim_type = currline(ind(end)+1:end-5);
    fclose(fid);


    stimTimes = A.data(:,2);
    out.stimTimes = stimTimes; % store the recorded timing of each stimulus record
    
    nValidStimFrames = nImages; %How many total frames actually imaged
    
    % Find average value of stimulus for each imaging frame.
    avrstimval = zeros(nValidStimFrames,1); 
    fstimval = zeros(nValidStimFrames,1);
    fstimpos1 = zeros(nValidStimFrames,1);
    fstimpos2 = zeros(nValidStimFrames,1);

    
    firstEntry = A.data(1,8);
     
    for k = firstEntry:nValidStimFrames
        inds = find(A.data(:,8) == k);
        if(~isempty(inds))
            stimval = A.data(inds,4);

            stimpos1 = A.data(inds,5); %ms, for jl type stimulus output files
            stimpos2 = A.data(inds,7); %ms, for jl type stimulus output files
            avrstimval(k) = mean(stimval);
            fstimval(k) = stimval(1);
            fstimpos1(k) = stimpos1(1);
            fstimpos2(k) = stimpos2(1);
            last_k_withStimEntries = k;
            
        %if scanning is faster than stimulus, use stimulus info of previous
        %frame that was written
        elseif(isempty(inds)) 
            inds = find(A.data(:,8) == last_k_withStimEntries); 
            stimval = A.data(inds,4);

            stimpos1 = A.data(inds,5); 
            stimpos2 = A.data(inds,7); 
            avrstimval(k) = mean(stimval);
            fstimval(k) = stimval(1);
            fstimpos1(k) = stimpos1(1);
            fstimpos2(k) = stimpos2(1);
        end
    end
    
    for k = 1:firstEntry-1
            avrstimval(k) = avrstimval(firstEntry);
            fstimval(k) = fstimval(firstEntry);
            fstimpos1(k) = fstimpos1(firstEntry);
            fstimpos2(k) = fstimpos2(firstEntry);
    end 

    save('stim', 'ch3', 'avrstimval', 'fstimval', 'frame_nums', ...
         'stim_type', 'fstimpos1', 'fstimpos2', 'stimTimes');

    out.ch3 = ch3;
    out.avrstimval = avrstimval;
    out.fstimval = fstimval;
    out.frame_nums = frame_nums;
    out.stim_type = stim_type;
    out.fstimpos1 = fstimpos1;
    out.fstimpos2 = fstimpos2;

%     save('data_file.mat', 'out');
    save('data_file.mat', 'out', '-v7.3'); 
    %% Make movie for Bruker images.
    if wantMovie
        % Raw data movie.
        makeTimeSeriesMovie(imageArray, 'original_series', 'mp4');
        % Aligned unfiltered data movie.
        makeTimeSeriesMovie(registeredImages, ...
                            'aligned_unfiltered_series_16bit', 'mp4');
        % Aligned median-filtered data movie.
        makeTimeSeriesMovie(medianFilteredImages, ...
                            'aligned_medianFiltered_series', 'mp4');
    end
end