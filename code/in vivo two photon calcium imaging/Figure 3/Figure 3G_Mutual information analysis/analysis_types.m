function neuronDataOfInterest = analysis_types(neuronData, options)
% Picks the analysus depending on the options

switch options.StimType
   
    case 'MI'
        switch options.MItype 
            case 'I_time_analysis'
                neuronDataOfInterest = MI_time_analysis(neuronData);
            case 'I_whole_trace' 
                neuronDataOfInterest = MI_calculation_randLum(neuronData);
        end
        
    otherwise
        neuronDataOfInterest = neuronData;
end

