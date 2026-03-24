

function out = aggregate_05StepsON_OFF(in)
% This function groups the different interesting epoch combinations in a
% Matrix.
% The epochs are distributed as:
% 
% 0 : 100% OFF
% 1 : 50% OFF
% 2 : GREY
% 3 : 50% ON
% 4 : 100% ON

% all interesting epoch-combinations (18 out of 25 possible ones)
N_EPOCHS = 20;

fprintf(' \n Remember: \n 0 : 100%% OFF \n 1 : 50%% OFF \n 2 : GREY \n 3 : 50%% ON \n 4 : 100%% ON \n\n');

%fi?gure out epochlength
epochlengths = zeros(1,length(in));

for ii=1:length(in)
    s = in(ii).istim;
    
    % Find beginning indices
    ds = diff(s);
    b_inds = find(ds)+1;
    duration = zeros(1,length(b_inds));
    
    % Calculate each epochduration for this ROI
    for kk =1:length(b_inds)-1
        duration(kk) = b_inds(kk+1)-b_inds(kk);
    end

    epochlengths(ii) = mode(duration);
end

epochlength = round(mean(epochlengths));
rats = zeros(length(in),N_EPOCHS,2*epochlength);

%% For each ROI
for ii=1:length(in)

     % df/f
     iratio = in(ii).iratio;
     is_100ON = (in(ii).istim == 4);
     begin = diff(is_100ON);
     b_inds100 = find(begin==1)+1;
     e_inds100 = b_inds100+epochlength-1;
     
     ON_100 = zeros(length(b_inds100),epochlength);
     
     for kk=1:length(b_inds100)
         % the last epochs might not be displayed fully, then jump out
         if(e_inds100(kk) > length(is_100ON))
             continue
         end
         
         ON_100(kk,:) = iratio(b_inds100(kk):e_inds100(kk));
     end
     
     ON_100(ON_100 == 0) = NaN;
     ON_100 = nanmean(ON_100,1);
     baseline = nanmean(ON_100(round(0.8*end):end));
     
     d_iratio = iratio./ baseline -1;

    
    % Now sort those epoch ratios in matrix "rats", for each epoch-combination 
    s = in(ii).istim;
    ds = diff(s);
     
    
    % Find indices for each trace
    % Problem: some combinations are left over (the ones which happen AFTER
    % a double-epoch was presented)
    b_inds = find(ds)+1;
    e_inds = b_inds+2*epochlength-1; % -1 because we want to go 'till the END
                                     % of the 2nd epoch.


    % Assign an identifying number to each epoch-combination in this ROI 
    % this is just an easier way to identify a combination using switch()
    for kk=1:length(b_inds) 
       % the last epochs might not be displayed fully, then jump out
       if(e_inds(kk) > length(s))
           continue
       end

       pairs(kk) = cantor(s(b_inds(kk)),s(e_inds(kk)));
    end
    
    % Store all examined pairs of this ROI in temp (approx 50)
    num_pairs = round(length(in(1).istim)/epochlength*2);
    temp = zeros(num_pairs,N_EPOCHS,2*epochlength);
    
    for jj=1:length(pairs)
        
        % If end of data is reached, jump out of loop.
        if(jj > numel(e_inds) || e_inds(jj) > length(d_iratio))
            continue;
        end
        
        b_ind = b_inds(jj);
        e_ind = e_inds(jj);
         
        % Check for interesting combinations. Combinations like (1,1) or (0,4)
        % are not interesting.
        switch pairs(jj)
            case cantor(0,1)% 100% OFF --> 50% OFF
                temp(jj,1,:)=d_iratio(b_ind:e_ind)';
            case cantor(0,2)% 100% OFF --> GREY
                temp(jj,2,:)=d_iratio(b_ind:e_ind)';
            case cantor(0,3)% 100% OFF --> 50% ON
                temp(jj,3,:)=d_iratio(b_ind:e_ind)';
            case cantor(0,4)% 100% OFF --> 100% ON
                temp(jj,4,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,0)% 50% OFF --> 100% OFF
                temp(jj,5,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,2)% 50% OFF --> GREY
                temp(jj,6,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,3)% 50% OFF --> 50% ON
                temp(jj,7,:)=d_iratio(b_ind:e_ind)';
            case cantor(1,4)% 50% OFF --> 100% ON
                temp(jj,8,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,0)% GREY --> 100% OFF 
                temp(jj,9,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,1)% GREY --> 50% OFF
                temp(jj,10,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,3)% GREY --> 50% ON
                temp(jj,11,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,4)% GREY --> 100% ON
                temp(jj,12,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,0)% 50% ON --> 100% OFF
                temp(jj,13,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,1)% 50% ON --> 50% OFF
                temp(jj,14,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,2)% 50% ON --> GREY
                temp(jj,15,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,4)% 50% ON --> 100% ON
                temp(jj,16,:)=d_iratio(b_ind:e_ind)';
            case cantor(4,0)% 100% ON --> 100% OFF
                temp(jj,17,:)=d_iratio(b_ind:e_ind)';   
            case cantor(4,1)% 100% ON --> 50% OFF
                temp(jj,18,:)=d_iratio(b_ind:e_ind)';    
            case cantor(4,2)% 100% ON --> GREY
                temp(jj,19,:)=d_iratio(b_ind:e_ind)';
            case cantor(4,3)% 100% ON --> 50% ON
                temp(jj,20,:)=d_iratio(b_ind:e_ind)';
            otherwise
                fprintf('No matching found for epoch combination (%d %d). Not interested? \n', cantorX(pairs(jj)),cantorY(pairs(jj)));
        end
    end
    
    temp(temp == 0) = NaN;
    rats(ii,:,:) = mean(temp,1,'omitnan');
    name{ii}=in(ii).name;
end   

% No baseline substraction
out.rats = rats;
out.neuron = [in.cell];
out.flyID = [in.flyID];
out.name = name;
