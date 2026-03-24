function imageArray = readStimTimeSeries(filename, imagingInfo)
    imagesFolder = regexp(filename, '.xml', 'split');
    
    cd(imagesFolder{1});
    % Initialize reader to the first stimulus frame.
    pathTo1stStimFile = which(imagingInfo.stimulusFrames{1}.fileName); 
    % Construct a Bio-Formats reader decorated with the Memoizer wrapper
    r = loci.formats.Memoizer(bfGetReader(), 0);
    % Initialize the reader with an input file to cache the reader
    r.setId(pathTo1stStimFile);
    % Close reader
    r.close()
    % If the reader has been cached in the call above, re-initializing the
    % reader will use the memo file and complete much faster especially for
    % large data
    r.setId(pathTo1stStimFile);
    % Perform additional work.
    % First file of the sequence is enough to read the whole time
    % series. Thus, we create a reader for the first file.
    reader = bfGetReader(pathTo1stStimFile);
    % Since the reader reads all images in the folder, set the
    % series to the stimulus series, i.e., iSeries = 2.
    iSeries = 2;
    reader.setSeries(iSeries - 1); % We need 0-based index.
    % Extract all the stimulus frames to a numeric array of dimensions
    % XYT.
    nImages = reader.getImageCount;
    imageArray = nan(reader.getSizeY, reader.getSizeX, reader.getSizeT);
    for jImage = 1 : nImages
        % Store all the nPlanes in the image array imageArray, it will
        % have dimensions nPixelsX x nPixelsY x nImages.
        imageArray(:, :, jImage) = bfGetPlane(reader, jImage);
    end
    % Close the readers
    reader.close()
    r.close()
end