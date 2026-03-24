function out = ROI_gcamp_analysis(in,roi_mode)

% Automatically choose ROIs

nframes = in.xml.frames;
fps = in.xml.framerate;

if(isfield(in,'AV1'))
    AV = in.AV1;
else
    AV = squeeze(sum(in.ch1a,3))/nframes; % The average image for ch1
end

if(nargin==1)||((nargin==2)&&(~roi_mode))
    
    figure;imagesc(AV);colormap gray;
    set(gca,'dataa',[1 1 1]);
    AV2=conv2(AV,ones(20)/400,'same');
    figure; imagesc(AV-AV2); colorbar;
    D = AV-AV2;

    DONE = 0;
    thresh = 0;
    but=1;
    md = imdilate(D > thresh,ones(3));
    mdr = imerode(md,ones(3));
    figure;
    while ~DONE
        threshkeep = thresh;
        buttonkeep = but;
        md = imdilate(D>thresh,ones(3));
        mdr = imerode(md,ones(3));

        subplot(2,1,1);
        hist(reshape(D,1,[]),200);
        set(gca,'yscale','log');
        subplot(2,1,2);
        imagesc(mdr);
        subplot(2,1,1);
        [thresh,dum,but] = ginput(1);
        DONE = (length(thresh)==0); % exit loop when return pressed...
    end
    mdrl = bwlabel(mdr);

    for i=1:max(max(mdrl))
        [ty,tx]=find(mdrl == i);
        cx(i)=mean(tx);
        cy(i)=mean(ty);
    end

    figure; imagesc(mdrl); hold on;
    plot(cx,cy,'ko');

    [x,y]=ginput;

    for i=1:length(x)
        d=(x(i)-cx).^2 + (y(i)-cy).^2;
        [dum,ch(i)]=min(d);
    end
    curLen = length(ch);
    
    % Repeat if necessary
    if (buttonkeep>1)
        figure;
        DONE = false;
        while ~DONE
            if(~isempty(thresh))
                threshkeep2 = thresh;                
            else
                thresh = threshkeep; 
            end
            md = imdilate(D>thresh,ones(3));
            mdr = imerode(md,ones(3));

            subplot(2,1,1);
            hist(reshape(D,1,[]),200);
            set(gca,'yscale','log');
            subplot(2,1,2);
            imagesc(mdr);
            subplot(2,1,1);
            [thresh,dum] = ginput(1);
            DONE = (length(thresh)==0); % exit loop when return pressed...
        end
        mdrl2 = bwlabel(mdr);

        for i=1:max(max(mdrl2))
            [ty,tx]=find(mdrl2 == i);
            cx(i)=mean(tx);
            cy(i)=mean(ty);
        end

        figure; imagesc(mdrl2); hold on;
        plot(cx,cy,'ko');

        [x,y]=ginput;

        for i=1:length(x)
            d=(x(i)-cx).^2 + (y(i)-cy).^2;
            [dum,ch(i+curLen)]=min(d);
        end
    end

    % generate masks
    nMasks = length(ch);
    masks = cell(length(ch),1);
    for i=1:length(ch)
        if(i<=curLen)
            curMask=imdilate(imdilate(mdrl == ch(i),ones(3)),ones(3)); % dilate
        else
            curMask=imdilate(imdilate(mdrl2 == ch(i),ones(3)),ones(3)); % dilate
        end
        masks{i} = curMask;
    end
    
elseif ((nargin==2)&&(~isnumeric(roi_mode)))

    figure;imagesc(AV);colormap gray;
    title('press enter when done selecting ROIs');
    done = 0;
    index = 1;
    while (~done)
        masks{index} = roipoly;
        done = waitforbuttonpress;
        index = index+1;
    end
    nMasks = index-1;

    curDir = pwd;
    cd(in.fileloc); % save masks in the T-series directory
    d = dir('curMasks*.mat');
    ind = length(d)+1;
    if (ind == 1)
        figure;imagesc(AV);colormap gray;
        title('select background region');
        NMask = roipoly;
    else
        load('curMasks1.mat','NMask');
    end
else
    curDir = pwd;
    cd(in.fileloc);
    load(sprintf('curMasks%d.mat',roi_mode));
    cd(curDir);
end

curDir = pwd;
cd(in.fileloc); % save masks in the T-series directory
if ((nargin==1)||(~isnumeric(roi_mode)))
    d = dir('curMasks*.mat');
    ind = length(d)+1;
    if (ind == 1)
        figure;imagesc(AV);colormap gray;
        title('select background region');
        NMask = roipoly;
    else
        load('curMasks1.mat','NMask');
    end
    save(sprintf('curMasks%d',ind),'masks','NMask','nMasks');
    disp(sprintf('saved curMasks%d',ind));
end
cd(curDir);

%% Generate ratio signals from all regions of interest - aligned data
   
out = in;
out.masks = masks;
out.NMask = NMask;
out.avSignal1 = zeros(nMasks,nframes);
out.dSignal1 = zeros(nMasks,nframes);
out.ratio = zeros(nMasks,nframes);
out.dRatio = zeros(nMasks,nframes);

if(~isfield(in,'AV1'))
    AV1 = squeeze(sum(in.ch1a,3))/nframes; % The average image for ch1
    out.AV1 = AV1;
end
if exist('in.xml.linesperframe') && exist('in.xml.linesperframe')
    masked = zeros(in.xml.linesperframe,in.xml.pixperline);
else
    masked = zeros(str2double(in.xml.linesPerFrame),str2double(in.xml.pixelsPerLine));
end
smask = zeros(nMasks,1);
for k = 1:nMasks
    smask(k) = sum(sum(masks{k}));
    masksi{k}= find(masks{k});
end
sNmask = sum(sum(NMask));
Nmaski = find(NMask);

for ind = 1:nframes
    A = double(squeeze(in.ch1a(:,:,ind)));
    for k = 1:nMasks

        masked = A(masksi{k});
        Nmasked = A(Nmaski);
        out.avSignal1(k,ind) = sum(masked)./smask(k);% summed signal in a ROI, normalized by ROI size
        out.dSignal1(k,ind) = out.avSignal1(k,ind) - sum((Nmasked))./sNmask; % background subtraction (by signal in background normalized by background ROI size)
    end
end

for i = 1:nMasks
    out.ratio(i,:) = out.avSignal1(i,:);
    out.dRatio(i,:) = out.dSignal1(i,:);
end

%% Bleaching analysis
% Analysis based on average signals, without background substraction

maxBleach = ones(2,1);
fdt = zeros(nMasks,1);
tdt = zeros(nMasks,1);

disp('bleaching info');
    disp(sprintf('neuron \t size \t ch1 b-rate \t ch2 b-rate \t ratio b-rate \t frac. dwell \t tot. dwell'));
for m = 1:nMasks
    % For every neuron
    bch = zeros(2,1);
    for ch = 1
        % For every channel
        avSig = eval(sprintf('out.avSignal%d(m,:)',ch));
        LavSig = log(avSig);
        p = polyfit(1:nframes,LavSig,1); 
        bch(ch) = p(1);

        if (p(1)<maxBleach(ch))
            maxBleach(ch) = p(1);
            nBleach = m;
        end
    end
    % For the ratio
    avSig = out.ratio(m,:);
    LavSig = log(avSig);
    p = polyfit(1:nframes,LavSig,1);
end     


%% Noise analysis

% filter
[b,a] = butter(2,0.1);

disp('noise info - with background substraction');
disp(sprintf(['neuron \t ch1 noise \t ch2 noise \t ratio noise \t\t signal \t SNR']));
for m = 1:nMasks
    % For every neuron
    noise = zeros(2,1);
    for ch = 1
        % For every channel
        dSig = eval(sprintf('out.dSignal%d(m,:)',ch));
        sdSig = filtfilt(b,a,dSig);
        ndSig = dSig-sdSig;
        noise(1,ch) = mean(ndSig);
        noise(2,ch) = std(ndSig);
    end
    % For ratio
    dSig = out.dRatio(m,:);
    sdSig = filtfilt(b,a,dSig);
    ndSig = dSig-sdSig;
    disp(sprintf('%d \t %0.2g+/-%0.3g \t %0.2g+/-%0.3g \t %0.3g \t %0.3g',m,noise(1,1),noise(2,1),mean(ndSig),std(ndSig),std(sdSig),std(sdSig)/std(ndSig))); %ms, only one channel for GCamP
end

disp('noise info - without background substraction');
disp(sprintf(['neuron \t ch1 noise \t ch2 noise \t ratio noise \t\t signal \t SNR']));
for m = 1:nMasks
    % For every neuron
    noise = zeros(2,1);
    for ch = 1
        % For every channel
        dSig = eval(sprintf('out.avSignal%d(m,:)',ch));
        sdSig = filtfilt(b,a,dSig);
        ndSig = dSig-sdSig;
        noise(1,ch) = mean(ndSig);
        noise(2,ch) = std(ndSig);
    end
    % For ratio
    dSig = out.ratio(m,:);
    sdSig = filtfilt(b,a,dSig);
    ndSig = dSig-sdSig;
    disp(sprintf('%d \t %0.2g+/-%0.3g \t %0.2g+/-%0.3g \t %0.3g \t %0.3g',m,noise(1,1),noise(2,1),mean(ndSig),std(ndSig),std(sdSig),std(sdSig)/std(ndSig)));
end

cd(curDir)
