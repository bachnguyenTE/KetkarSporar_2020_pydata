% ********************************************************************
% ********************************************************************
% ********************************************************************
%        Random luminance steps mutual information analysis
% ********************************************************************
% ********************************************************************
% ********************************************************************

% please add function path on line number 16
%% Initialize
clc;
clear all;
close all;

% Add function path
addpath(genpath('/your_functionPath'))


% Pre-defined parameters
GenotypeNumber = 2;
options.interpolationRate = 10;
options.StimType = 'MI'; % Determines the analysis in the upcoming functions
options.MItype = 'I_whole_trace' ; % 'I_whole_trace' | 'I_peak_late' | 'I_time_analysis'
options.Database = 'Katja' % Using Katja Sporar's database
options.NegCorr = 1; % Using negative correlation thresholding of ROIs
% open indexStorage % Here are the boolean operations for selecting indices

%% Data acquisiton and pre-processing

for iGenotype = 1:GenotypeNumber
    % Acquire & pre-process data
    neuronDataOfInterest{iGenotype}  = data_acquire_pre_process ...
    ( options, 1) ;

    % Give a name to the current criteria for plotting
    genotypeName{iGenotype} = ...
        inputdlg('Enter the name of the criteria');
end


switch options.MItype
    case 'I_whole_trace'
        [statsStruct] = ...
            MI_whole_trace_analysis_plot(GenotypeNumber, neuronDataOfInterest, ...
            genotypeName);
    case 'I_time_analysis'
        MI_time_analysis_plot(GenotypeNumber, neuronDataOfInterest, genotypeName);
end


%% Saving data
answer = questdlg('Do you want to save the data?');
if strcmp(answer,'Yes')
    for iGenotype = 1:GenotypeNumber

        % -- Adding cell type data -- 
        fieldNames = fieldnames(neuronDataOfInterest{iGenotype});

        % Size of the data
        currSize = max(size(neuronDataOfInterest{iGenotype}.(fieldNames{1})));

        % Create a cell for cell type data
        currCellType = cell(currSize,1);
        for iFill = 1:currSize
            currCellType{iFill} = char(genotypeName{iGenotype});
        end

        % Add cell type data
        neuronDataOfInterest{iGenotype}.CellType  = currCellType;

         % -- Saving data -- 
         struct2csv(neuronDataOfInterest{iGenotype})


    end
end






