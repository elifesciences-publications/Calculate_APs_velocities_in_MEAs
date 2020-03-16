function  my_electrode_indices = Aux_Electrode_Index_From_Label( my_electrode_labels, electrode_labels )
% Get electrodes indices from electrodes labels
%
%   Syntax:
%       Aux_Electrode_Index_From_Label( my_electrode_labels )
%
%   Input variables:
%       my_electrode_labels: list of labels, e.g {'O9'; 'O10'; 'O11'; 'O12'} 
%
%   Output parameters:
%       my_electrode_indices: 1D array with selected electrodes ids
%
%   Requires:
%       electrode_labels 
%
%   Provides:
%       n.a.
%
%   Example:
%       my_electrode_indices = Aux_Electrode_Index_From_Label( {'O9'; 'O10'} )
%
% Jose Mateus, Paulo Aguiar
% INEB/i3S, Mar 2018
% pauloaguiar@ineb.up.pt
% -----------------------------------------------------------------------
   
    if ~iscell( my_electrode_labels )
        my_electrode_labels = { my_electrode_labels };
    end

    N = numel( my_electrode_labels );
    my_electrode_indices = zeros( N, 1 );    
    
    for i = 1:N
        temp = find( strcmp( electrode_labels, my_electrode_labels{i} ) ); % ensure correct order of electrode IDs (indexes) 
        if isempty( temp )
            disp(['ERROR: non valid electrode label - ', my_electrode_labels{i} ]);
            my_electrode_indices = -1;
            beep
            return;
        else
            my_electrode_indices(i) = temp;
        end
    end    

end