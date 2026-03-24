function out = aggregate_standingStripe_ks_grey_baseline(in)


% figure out how many epochs there are
s = in(1).istim;
n_epoch = sum(unique(s)>0);  %count number of unique epoch numbers but ignore 0


epochlength = zeros(length(in),1);
zeroInterleaveLength = zeros(length(in),1);
Combined_epochlength = zeros(length(in),1);
% loop over all of the stimuli to find the epochlength to be used:
for ii=1:length(in)
    s = in(ii).istim;
    % finding the actual epochs
    ds = diff(s); 

    % begining of the epoch, ie when the stimulus steps back down to 0
    % (which is the interstimulus blank period)
    binds = find(ds<0)+1; % ignoring the first epoch, finding the beginning of the inter-epoch interval (0s in stim strace) 
    einds = find(ds>0)+1; % finding the beginning of the actual epochs
    
    % find out which is shorter so that subtraction works out
    len = min(length(binds),length(einds));
    
    % duration of the numbered epoch
    dinds = binds(1:len)-einds(1:len);
    %ms: automatically determine epoch duration
    epochlength(ii)=round(mean(dinds));
    
    % duration of the zero interleave epochs
    zinds = einds(2:len)- binds(1:len - 1);
    zeroInterleaveLength(ii) = round(mean(zinds));
    
    %automatically determine epoch duration, by looking for # of frames
    %between each step from a number (ie 8#) to 0. so includes both 000000
    %and 8888888 period for example
    Combined_epochlength(ii)=round(mean(diff(binds)));
    
end

% duration of the numbered epochs '555555555' or '22222222' 
epoch_dur = round(mean(epochlength));% if frame rate was ~10Hz and using 1 sec StandingStripe this value should be around 10.
% This duration of the zero epochs '00000000' 
zero_dur = round(mean(zeroInterleaveLength)); % if frame rate was ~10Hz and using 1 sec StandingStripe this value should be around 10.

% full duration include preceding zero epoch, 'real' epoch and the follwing
% zero epoch
full_dur = zero_dur + epoch_dur + zero_dur;
stimstruct = zeros(full_dur,1);
stimstruct(zero_dur+1:zero_dur+epoch_dur) = 1; %filling in the ones


rats = zeros([length(in),n_epoch,full_dur]);
name = cell(length(in),1);


% extract the data as organized by epochs and store in for further processing
for ii=1:length(in)
     % df/f
     iratio = in(ii).iratio;
     is_grey = (in(ii).istim == 0);
     begin = diff(is_grey);
     b_indsgrey = find(begin==1)+1;
     e_indsgrey = b_indsgrey+zero_dur-1;
     
     
     grey = zeros(length(b_indsgrey),zero_dur);
     
     for kk=1:length(b_indsgrey)
         % the last epochs might not be displayed fully, then jump out
         if(e_indsgrey(kk) > length(is_grey))
             continue
         end
         
         grey(kk,:) = iratio(b_indsgrey(kk):e_indsgrey(kk));
     end
     
     grey(grey == 0) = NaN;
     grey = nanmean(grey,1);
     baseline = nanmean(grey(round(0.5*end):end));
     
     d_iratio = iratio./ baseline -1;
     

    s = in(ii).istim;
    % finding the actual epochs
    ds = diff(s);
    inds = find(ds>0)+1; % 
    inds = inds(2:end-1); % clipping off the first and the last epoch
    
    epochs = s(inds);
    
    pos = in(ii).ifstimpos1;
    
    for jj=1:n_epoch % number of epochs
        
        I = find(epochs==jj); %ms: take all the epochs
        bind = inds(I);
        
        %make variable to store the traces
        allEpochs = zeros(length(bind),full_dur);
        for mm = 1:length(bind)
            i_beginingOfTrace = bind(mm)- zero_dur;
            i_endOfTrace = bind(mm) + epoch_dur + zero_dur -1;
            % extract full part of the trace that we want, includes prior zero
            % period, the actuall epoch and following zero period
            allEpochs(mm,:) = d_iratio(i_beginingOfTrace:i_endOfTrace);

            temp = mean(allEpochs,1);
        end
        
        temp = temp- mean(temp(1:zero_dur));
  
        rats(ii,jj,:) = temp;
    end
    name{ii}=in(ii).name;

end

% save values into out variable and return to function.
out.rats = rats;
out.neuron = [in.cell];
out.flyID = [in.flyID];
out.name = name;
out.stimstruct = stimstruct;

% end

     
     
     
     
     
     
