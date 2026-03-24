%% This script predicts contrast-sensitive LMC responses for tested OFF edges
% based on the data in Laughlin et al 1987, with the assumption that the fly 
% adapts to mean luminance, and plots these predictions against background luminance.

clearvars
%% fit of the contrast-response curve from Laughlin et al, 1987
data=csvread('contrastResponses.csv'); % response per unit contrast
f=fit(data(:,1),data(:,2),'a/(1+exp(-b*(x-c)))');
figure; plot(f,data(:,1),data(:,2));

%% define stimulus parameters
Ilevels=0:15; %LED levels
fractionLevels=Ilevels/15;
ndVals=[0 0.9 1.8 2.7 3.6]; %neutral density filter values
ndNames={'noND','ND0.9','ND1.8','ND2.7','ND3.6'};
maxL=1176758.88; %max native LED luminance in photons per receptor per s.
nativeLs=[maxL/15,maxL*2/15,maxL*4/15,maxL*8/15,maxL]; 
          %corresponding to the LED levels 1, 2, 4, 8, 15.
lumValues=[nativeLs;... %for no ND
    nativeLs./8;... %for ND 0.9
    nativeLs./64;... %for ND 1.8
    nativeLs./512;... %for ND 2.7
    nativeLs./4096]; %for ND 3.6
ndMeanAdaptLevels=(lumValues(:,5)/15)*3; 
                  %corresponding to LED level 3 (~stimulus mean) in each ND.
bglevel=fractionLevels(1+[1,2,4,8,15]); %background LED levels
offlevel=[0,0,0,0,0]; %OFF edge
contsbg=(bglevel-((bglevel(5)/15)*3))/((bglevel(5)/15)*3);
              %contrast between background lum and adapting lum
contsOff=(offlevel-((bglevel(5)/15)*3))/((bglevel(5)/15)*3);
              %contrast wrt adapting lum

%% calculate contrast-sensitive LMC responses to OFF edges
contrastResponse=f.a./(1+exp(-f.b*(log10(ndMeanAdaptLevels)-f.c)));
contrastResponse=contrastResponse/100; %percentage
% Now, calculate responses to the contrast between preStim and adapting luminance,
% as well as postStim and adapting luminance. These responses lie on the  
% contrast vs response curve centered at 0 contrast where 'contrastResponse' is slope
LMCbg=1./(1+exp(contrastResponse*(contsbg))); %negative sign canlceled by the negative slope for L2
LMCoff=1./(1+exp(contrastResponse*(contsOff))); %negative sign canlceled by the negative slope for L2

LMCresponse=LMCoff-LMCbg; %response to OFF edges
% each row of the matrices in this section: one nd filter