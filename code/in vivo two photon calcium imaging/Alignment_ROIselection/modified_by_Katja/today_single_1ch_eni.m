clc;
clear all;
close all;

% Single file analysis


i = 4;
file = ['/Users/ksporar/2p/191202ks_fly4/Image ' num2str(i)];

disp(file);

curDir = pwd;cd(file)

m = dir('data_file*');
    
if(~isempty(m))

    load('data_file');

    out = ROI_gcamp_analysis(out, true); %true if you want to select ROIs
    out = FFFlash_res_display_1ch(out, 2);

    save_processed_data_1ch_eni(out);

    d = dir('_stimulus_*');
    fid = fopen(d.name,'r');
    currline = fgetl(fid);
    ind = strfind(currline,'\');
    disp(sprintf('stimulus: %s',currline(ind(end)+1:end)));
    fclose(fid);

else
    m = dir('*pData.mat');
    load(m.name);
    out = FFFlash_res_display_1ch(strct,2);
    ind = strfind(out.fileloc,'/');
    ind2 = strfind(file,'/');
    out.fileloc = [file(1:ind2(end)) out.fileloc((ind(end)+1):end)];
    save_processed_data_1ch(out);
    
end