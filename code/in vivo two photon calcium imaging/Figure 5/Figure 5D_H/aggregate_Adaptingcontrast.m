

function out = aggregate_Adaptingcontrast(in,OFF_stimulus)
% This function groups the different interesting epoch combinations in a
% Matrix.

N_EPOCHS = 7;
N_APOCHS = 1;
fprintf(' \n Remember: \n 0 : 100%% OFF \n 1 : 50%% OFF \n 2 : GREY \n 3 : 50%% ON \n 4 : 100%% ON \n\n');

%figure out epochlength
epochlengths1 = zeros;
epochlengths2 = zeros;
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
    %duration = sort(duration);
    for kk = 1:(length(b_inds)/2)
        if kk == 1
            epochlengths1(ii) = duration (2);
            epochlengths2(ii) = duration (1);
            
        else
            
    epochlengths1(kk) = (duration(2*kk));
    epochlengths2(kk) = (duration(2*kk-1));
        end
    end
    
    
    


% For each ROI (length(in)) the average of

realepochlength1 = 300;
realepochlength2 = 60;
if ii == 1
rats = zeros(length(in),N_EPOCHS,(mean(realepochlength1)+mean(realepochlength2)));
end





    % df/f
    iratio = in(ii).iratio;
    istumulilog=in(ii).istim==0;
    lowest10 = iratio(istumulilog==1);
    lowest10 = sort (lowest10);
    
    lowest10 = lowest10(1,round(0.25*end):round(0.75*end));
    iratio = iratio/mean(lowest10)-1;
    d_iratio = iratio;
    
    
    % Find indices for each trace
    % Problem: some combinations are left over (the ones which happen AFTER
    % a double-epoch was presented)
    e_inds = zeros(1,length(b_inds));
    for kk = 1:length (b_inds);
    b_inds = find(ds)-149;
    e_inds(kk) = ((b_inds(kk)+epochlengths1(round(0.5*kk)))+(epochlengths2(round(0.5*kk)))-1); % -1 because we want to go 'till the END
        
    end
    for kk=1:length(b_inds) 
       % the last epochs might not be displayed fully, then jump out
       if((e_inds(kk)) > length(s))
           continue
       end

       pairs(kk) = cantor(s((b_inds(kk)+150)),s((e_inds(kk))));
    end
    
    
     
    % Store all examined pairs of this ROI in temp (approx 50)
    num_pairs = round(length(in(1).istim)/((realepochlength1)+(realepochlength2)));
    temp = zeros(num_pairs,N_EPOCHS,((realepochlength1)+(realepochlength2)));

    for jj=1:length(pairs)
        
        % If end of data is reached, jump out of loop.
        if(jj > numel(e_inds) || e_inds(jj) > length(d_iratio))
            continue;
        end
        
        b_ind = b_inds(jj);
        e_ind = e_inds(jj);
        if length(d_iratio(b_ind:e_ind))>360;
            e_ind = e_ind - 1;
        elseif length(d_iratio(b_ind:e_ind))<360;
            e_ind = e_ind + 1;
        end

        if OFF_stimulus == 1 && e_ind - b_ind ==359; 
            switch pairs(jj)
            case cantor(1,0)% 100% OFF --> 50% OFF
                temp(jj,1,:)=d_iratio(b_ind:e_ind)';
            case cantor(2,0)% 100% OFF --> GREY
                temp(jj,2,:)=d_iratio(b_ind:e_ind)';
            case cantor(3,0)% 100% OFF --> 50% ON
                temp(jj,3,:)=d_iratio(b_ind:e_ind)';
            case cantor(4,0)% 100% OFF --> 100% ON
                temp(jj,4,:)=d_iratio(b_ind:e_ind)';
            case cantor(5,0)% 50% OFF --> 100% OFF
                temp(jj,5,:)=d_iratio(b_ind:e_ind)';
            case cantor(6,0)% 50% OFF --> GREY
                temp(jj,6,:)=d_iratio(b_ind:e_ind)';
            case cantor(7,0)% 50% OFF --> 50% ON
                temp(jj,7,:)=d_iratio(b_ind:e_ind)';
            
            otherwise
                fprintf('No matching found for epoch combination (%d %d). Not interested? \n', cantorX(pairs(jj)),cantorY(pairs(jj)));
           end
        end
        
        
    
    temp(temp == 0) = NaN;
    rats(ii,:,:) = nanmean(temp,1);
    name{ii}=in(ii).name;
        
    end
   
end   

% No baseline substraction
out.rats = rats;
out.neuron = [in.cell];
out.flyID = [in.flyID];
out.name = name;
