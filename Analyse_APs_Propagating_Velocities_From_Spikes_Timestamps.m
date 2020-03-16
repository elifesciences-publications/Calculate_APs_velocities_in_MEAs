% Analysis script for calculating action potentials propagation velocities
% using spikes timestamps in specified electrodes in a MEA device (coupled
% with microfluidics)
% Jose Mateus, Miguel Aroso, Paulo Aguiar
% INEB/i3S, Mar 2019
% pauloaguiar@ineb.up.pt
 
%% Load data
load('SampleData_SpikesTimes_and_ElecLabels.mat');
 
%% Variables
elec_labels = 'ABCDEFGHJKLMNOPR';
elec_num    = 9:15;

%% Pre-allocate
% for e = 1:numel( elec_labels )
%     eval( ['prop_vel_antero_', elec_labels(e), ' = struct([]);'] );
%     eval( ['prop_vel_retro_' , elec_labels(e), ' = struct([]);'] );
% end

%% Save Propagating spikes per microchannel & Propagation Velocities
for e = 1:numel( elec_labels )
    
    disp( ['Detecting propagating spikes, and calculating their velocity, in microchannel ', elec_labels(e)] );

    % construct variable my_electrodes
    my_electrodes = cell( numel(elec_num), 1 );
    for n = 1:numel(elec_num)
        my_electrodes{n} = [ elec_labels(e), num2str( elec_num(n) ) ];
    end
    
    % detect propagation and calculate velocity
    [prop_spikes_antero, prop_spikes_retro, prop_vel_antero, prop_vel_retro] = Calculate_Propagating_Spikes( my_electrodes, 100, 1, 3, electrode_labels, spike_times_elec_ms);

    % organize results in Workspace
    eval( ['prop_spikes_antero_', elec_labels(e), '= prop_spikes_antero;'] );     
    eval( ['prop_vel_antero_',    elec_labels(e), '= prop_vel_antero;'] );
    eval( ['prop_spikes_retro_',  elec_labels(e), '= prop_spikes_retro;'] );
    eval( ['prop_vel_retro_',     elec_labels(e), '= prop_vel_retro;'] );
    
    % save data
    % if numel(prop_spikes_antero) > 0
    %     save( ['Microchannel_', elec_labels(e), '_antero_spikes.mat'], 'prop_spikes_antero');
    %     save( ['Microchannel_', elec_labels(e), '_antero_PV.mat'],     'prop_vel_antero');
    % end
    % 
    % if numel(prop_spikes_retro) > 0
    %     save( ['Microchannel_', elec_labels(e), '_retro_spikes.mat'], 'prop_spikes_retro');
    %     save( ['Microchannel_', elec_labels(e), '_retro_PV.mat'],     'prop_vel_retro');
    % end

end

clear ('prop_spikes_antero', 'prop_spikes_retro', 'prop_vel_antero', 'prop_vel_retro', 'e', 'n');
