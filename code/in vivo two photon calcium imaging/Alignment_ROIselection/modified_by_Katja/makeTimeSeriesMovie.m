function makeTimeSeriesMovie(timeSeries, movieFilename, movieFormat)
    if nargin < 3
        movieFormat = 'mp4';
    end
    switch movieFormat
        case or('mp4', 'm4v')
            profile = 'MPEG-4';
        case 'avi'
            profile = 'Grayscale AVI';
        case 'mj2'
            profile = 'Archival';
        otherwise
            movieFormat = 'mp4';
            profile = 'MPEG-4';
    end
            
    writerObj = VideoWriter([movieFilename '.' movieFormat], profile);
    open(writerObj);
    for iFrame = 1 : size(timeSeries, 3)
        imagesc(squeeze(timeSeries(:, :, iFrame)))
        frame = getframe;
        writeVideo(writerObj, frame);
    end
    close(writerObj);
    close(gcf)
end