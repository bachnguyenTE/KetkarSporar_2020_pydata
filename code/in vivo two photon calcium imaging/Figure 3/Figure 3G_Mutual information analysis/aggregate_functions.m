function neuronData = aggregate_functions(neuronRawData, options)
% Picks the aggregate depending on the data
% Deleted all the other cases and left just the relevant one for
% Ketkar_Sporar et al.
switch options.StimType
    
    case 'MI'
        neuronData = neuronRawData ;
 
end
