%% Perform alignment based on maximizing image crosscorrelation in the Fourier space.
function [registeredImages, unregisteredImages, ...
          medianFilteredImages] = fourierCrossCorrelAlignment(imageArray, ...
                                                              refFrame, ...
                                                              filetype, ...
                                                              savedFileName)
    if nargin < 3 || isempty(filetype)
        filetype = 'xml';
    end
    dimension = size(imageArray);
    height = dimension(1);
    width = dimension(2);
    nImages = dimension(3);
    % Initialize image arrays:
    registeredImages = zeros(height, width, nImages, 'uint16');
    unregisteredImages = zeros(height, width, nImages, 'uint16');
    medianFiltImages = zeros(height * width, nImages);
    % Fourier transform of reference image, DC in (1,1)   [DO NOT FFTSHIFT].
    refImageFT = fft2(refFrame);
    % Upsampling factor (integer). Images will be registered to within 1/usfac 
    % of a pixel. For example usfac = 20 means the images will be registered 
    % within 1/20 of a pixel. (default = 1)
    upSamplingFactor = 1;

    for jImage = 1: nImages
        % Fourier transform of image to register, DC in (1,1) [DO NOT FFTSHIFT].
        switch filetype
            case 'lsm'
                jImageToRegister = imageArray(jImage).data;
                jImageFT = fft2(im2double(jImageToRegister));
            case 'xml'
                jImageToRegister = imageArray(:, :, jImage);
                jImageFT = fft2(im2double(jImageToRegister));
            otherwise
                error(['Unrecognized file type: ' filetype])        
        end

        [~, registeredImageFT] = dftregistration(refImageFT, jImageFT, ...
                                                 upSamplingFactor);    
        % Convert image to 16 bit format. Check this for use with different bit
        % depths.
        registeredImages(:, :, jImage) = uint16(ifft2(registeredImageFT));
        unregisteredImages(:, :, jImage) = jImageToRegister;

        if jImage == 1
            imwrite(registeredImages(:, :, jImage), savedFileName);
        else
            imwrite(registeredImages(:, :, jImage), savedFileName, 'WriteMode', 'append');
        end

        % Median filter (needed for normalization). It was set to 5 pixels
        % neighbourhood.
        neighbourhoodSize = 3 * [1 1];
        jImageMedianFilter = medfilt2(registeredImages(:, :, jImage), ...
                                      neighbourhoodSize);
        medianFiltImages(:, jImage) = reshape(jImageMedianFilter, ...
                                              height * width, 1) + 1;
    end

    % Reshape median-filtered images into a time series of images.
    medianFilteredImages = reshape(medianFiltImages, ...
                                   [height, width, nImages]);
end