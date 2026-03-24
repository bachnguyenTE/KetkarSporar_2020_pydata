function out = FFFlash_res_display_1ch(in,dur)
%% Find stimulus timing - modify after clean-up script is written

out = in;
nframes = in.xml.frames;
fps = in.xml.framerate;

thresh = 0.5;
mask = out.fstimval>thresh;
dmask = mask(2:end)-mask(1:end-1);
dmask = [0; dmask];

%% Plot ratio signals from ROIs together with stimulus timing

nMasks = length(in.masks);
cm = colormap('lines');

for i = 1:nMasks
    plot((1:nframes)/fps,in.avSignal1(i,:),'color',cm(i,:));hold on;
end
xlabel('time (sec)');
title('Signal in region of interest - before background substraction, aligned data');

figure; 
for i = 1:nMasks
    plot((1:nframes)/fps,(in.dSignal1(i,:)/mean(in.dSignal1(i,:)))+1*(i-1),'color',max(cm(i,:),.7));hold on;
    plot((1:nframes)/fps,in.dRatio(i,:)/mean(in.dRatio(i,:))+1*(i-1),'color',cm(i,:),'linewidth',2);
end
plot((1:nframes)/fps,out.fstimval*0.2,'LineWidth',2);
axis([0 nframes/fps 0 i+1]);

inds = find(dmask~=0);
for k = 1:length(inds)
    if(dmask(inds(k))>0)
        line([inds(k)/fps inds(k)/fps],[0 nMasks+1],'color',cm(i,:),'LineWidth',1,'LineStyle','-');
    else
        line([inds(k)/fps inds(k)/fps],[0 nMasks+1],'color',cm(i,:),'LineWidth',1,'LineStyle','--');
    end
end
xlabel('time (sec)');
title('Ratio, YFP, CFP signals - aligned data');

% Create a colored map of ROIs
if exist('in.xml.linesperframe') && exist('in.xml.pixperline')
    CMask = zeros(in.xml.linesperframe, in.xml.pixperline, 3); % Luis 13.11.2015
else
    CMask = zeros(str2double(in.xml.linesPerFrame), str2double(in.xml.pixelsPerLine), 3);
end

for i = 1:nMasks
    curColor = cm(i,:);
    curMask = cat(3,curColor(1).*flipud(in.masks{i}),curColor(2).*flipud(in.masks{i}),curColor(3).*flipud(in.masks{i}));
    CMask = CMask + curMask;
end

if(~isfield(in,'AV1'))
    AV = squeeze(sum(in.ch1a,3))/nframes; % The average image
    AV = im2double(AV);
    AV = AV./max(AV(:));
else
    AV = in.AV1;
end
figure;imshow(flipud(AV),[]);
hold on;h = imshow(CMask);
set(h,'AlphaData',0.5);
