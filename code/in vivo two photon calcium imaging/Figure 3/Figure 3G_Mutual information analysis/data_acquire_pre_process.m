function neuronDataOfInterest  = data_acquire_pre_process...
    ( options, ASKCriteria )
% This function asks for experimental groups and processes the data
% Input: 
%   options : structure
%       A structure with analysis options as fields.
%   ASKCriteria : bool
%       If user wants to manually enter the criteria for genotype
%       selection. see indexStorage.m for examples.

%%%%%%%%%%%%%%%% ---- START FUNCTION ---- %%%%%%%%%%%%%%%%
% open indexStorage --> for copy pasting flies of interest

% if ASKCriteria are not given make use the pre-defined criteria
if nargin == 1
    ASKCriteria = 0;
end
% add the path of needed functions
addpath(genpath('/your_functionPath'))

% Look for processed database (left only Katja's database which is relevant
% for the paper)
switch options.Database 
    case 'Katja'
         % Locate the Master Folder
        pdatapath='/Volumes/HD-SP1/Katja_pData';
        addpath(pdatapath);
        try
            cd(pdatapath);
        catch error1
            error('Katja''s pData path not accesible')
        end

        % Run database 
        database_select_samples_ks;
        
end
        




%% Data acquisiton and pre-processing    
% Criteria selection via input
if ASKCriteria
    try
        criteria = input('Please enter the criteria you want \n ->');
    catch error2
        errordlg('Criteria are not correct')
        error2
        return
    end
    
 
    indicesOfInterest  = find(criteria); %For prompting the criteria

    fprintf('%d criteria are found \n', length(indicesOfInterest))
    if isempty(indicesOfInterest)
        errordlg('No criteria found please try again')
        return
    end

else
    errordlg('Automatic genotype selection is not yet available for the current version.')
end

%% Pre-processing 

% Interpolation depending on stimulus (just left paper methods)
switch options.Database 
    case 'Katja'
         % Locate the Master Folder
         
        neuronStructure = create_neuron_structure_all_ks(indicesOfInterest);
        neuronRawData = load_neuron_data10Hz_ks(neuronStructure,pdatapath);
        
end

%% Choosing the correct aggregate function 

neuronData = aggregate_functions(neuronRawData, options);
if strcmp(options.StimType,'MI') || strcmp(options.StimType,'ND_analysis_random_steps')
    if options.NegCorr
        neuronData(1).negCorrROIs = rand_steps_corr_katja_method(neuronRawData);
    end
end
%% Choosing the correct analysis for stimulus Type

neuronDataOfInterest = analysis_types(neuronData, options);

