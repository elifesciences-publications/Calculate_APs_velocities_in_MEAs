function [prop_spikes_antero_ms, prop_spikes_retro_ms, prop_vel_antero, prop_vel_retro] = Calculate_Propagating_Spikes( my_electrode_labels, elec_distance_um, max_inter_elec_delay_ms, min_inter_spike_delay_ms, electrode_labels, spike_times_elec_ms )
% Detects propagating spikes along an electrode sequence (my_electrode_labels) based on time-delays
% Plots a specialized raster plot with propagating spikes being color coded 
%
%   Syntax:
%       [prop_spikes_antero_ms, prop_spikes_retro_ms, prop_vel_antero, prop_vel_retro] = Calculate_Propagating_Spikes( my_electrode_labels, elec_distance_um, max_inter_elec_delay_ms, min_inter_spike_delay_ms )
%
%   Input parameters:
%       my_electrode_labels: electrode sequence to analyze (e.g {'O9'; 'O10'; 'O11'; 'O12'})
%       elec_distance_um: inter-electrode distance in um (e.g 100.0)
%       max_inter_elec_delay_ms: maximum time delay in ms between consecutive electrodes to consider a propagating spike
%       min_inter_spike_delay_ms: minimum time delay in ms between consecutive spikes, to avoid calculation problems due to bursts
% 
%   Output variables:
%       propagating_spikes_antero: array with timestamps of antero propagating spikes for the electrode sequence 
%       propagating_spikes_retro: array with timestamps of retro propagating spikes for the electrode sequence 
%
%   Requires:
%       obj.spike_times_elec_ms
%
%   Provides:
%       n.a.
%
% Jose Mateus, Ana Geros, Miguel Aroso, Paulo Aguiar
% INEB/i3S, Mar 2019
% pauloaguiar@ineb.up.pt
% -----------------------------------------------------------------------

    if numel(my_electrode_labels) < 2       % min 2 electrodes to detect propagating spikes
       disp('ERROR!: Not enough electrodes selected');
       return
    end

    disp({'WARNING:'; 'Selected electrodes are assumed to be presented in spatial order!'; '(with electrodes closest to soma presented first)'}); 

    
    % Get_Electrode_Index_From_Label 
    electrodes_indices = Aux_Electrode_Index_From_Label( my_electrode_labels, electrode_labels );
    N_elec = numel( my_electrode_labels );    
    
    
    % Prealocate array for selected spike times - use number of spikes in first electrode 
    prop_spikes_antero_ms = NaN * ones( N_elec, numel( spike_times_elec_ms{ electrodes_indices(1) } ) );
    prop_spikes_retro_ms = NaN * ones( N_elec, numel( spike_times_elec_ms{ electrodes_indices(1) } ) );
    
    
    % Use first elect as refernces, and fill array with closest spikes, in both directions
    for spk = 1 : numel( spike_times_elec_ms{ electrodes_indices(1) } )
        
        Tspike_ANT_ms = spike_times_elec_ms{ electrodes_indices(1) }(spk); 
        Tspike_RET_ms = Tspike_ANT_ms; 
        prop_spikes_antero_ms(1,spk) = Tspike_ANT_ms;
        prop_spikes_retro_ms(1,spk) = Tspike_RET_ms;
        
        for elec = 2 : numel(electrodes_indices)
            
            targets_ms = spike_times_elec_ms{ electrodes_indices(elec) };
            targets_ANT_ms = targets_ms( targets_ms > Tspike_ANT_ms );
            targets_RET_ms = targets_ms( targets_ms < Tspike_RET_ms );
            
            % take care of ANT
            diffs_ANT_ms = targets_ms - Tspike_ANT_ms;
            diffs_ANT_ms( diffs_ANT_ms < 0 ) = Inf;
            [delay_ms, ind] = min( diffs_ANT_ms );            
            if ~isempty( targets_ANT_ms ) && delay_ms < max_inter_elec_delay_ms
                Tspike_ANT_ms = targets_ANT_ms(1);
                % test for condition of being part of a burst
                if ind > 1
                    if Tspike_ANT_ms - targets_ms(ind-1) < min_inter_spike_delay_ms
                        Tspike_ANT_ms = NaN;
                    end
                end
                if ind < numel( targets_ms )
                    if targets_ms(ind+1) - Tspike_ANT_ms < min_inter_spike_delay_ms
                        Tspike_ANT_ms = NaN;
                    end
                end 
                prop_spikes_antero_ms(elec,spk) = Tspike_ANT_ms;
            end
            
            % take care of RET
            diffs_RET_ms = targets_ms - Tspike_RET_ms;
            diffs_RET_ms( diffs_RET_ms > 0 ) = -Inf;
            [delay_ms, ind] = max( diffs_RET_ms );
            if ~isempty( targets_RET_ms ) && -delay_ms < max_inter_elec_delay_ms
                Tspike_RET_ms = targets_RET_ms(end);
                % test for condition of being part of a burst
                if ind > 1
                    if Tspike_RET_ms - targets_ms(ind-1) < min_inter_spike_delay_ms
                        Tspike_RET_ms = NaN;
                    end
                end
                if ind < numel( targets_ms )
                    if targets_ms(ind+1) - Tspike_RET_ms < min_inter_spike_delay_ms
                        Tspike_RET_ms = NaN;
                    end
                end 
                prop_spikes_retro_ms(elec,spk) = Tspike_RET_ms;
            end
            
        end
        
        % flag all non-ordered sequences
        if ~issorted( prop_spikes_antero_ms(:,spk), 'ascend'  )
            prop_spikes_antero_ms(1,spk) = NaN;
        end
        if ~issorted( prop_spikes_retro_ms(:,spk), 'descend' )
            prop_spikes_retro_ms(1,spk) = NaN;
        end
        
    end
    
    
    %% Remove all sequences which are not ordered
    bad_seqs = isnan( prod( prop_spikes_antero_ms, 1 ) ); % workaround since issorted([1, nan]) is true
    prop_spikes_antero_ms(:,bad_seqs) = [];    
    bad_seqs = isnan( prod( prop_spikes_retro_ms, 1 ) ); % workaround since issorted([1, nan]) is true
    prop_spikes_retro_ms(:,bad_seqs) = [];    

       
    %% Calculation of propagation velocities by two methods:
    % 1) through dx/dt given the extreme values
    total_dist_um = elec_distance_um * ( N_elec - 1 );
    prop_vel_antero.first_last_elec_m_per_sec = 1.0e-3 * total_dist_um ./ ( prop_spikes_antero_ms(N_elec,:) - prop_spikes_antero_ms(1,:) );
    prop_vel_retro.first_last_elec_m_per_sec  = 1.0e-3 * total_dist_um ./ ( prop_spikes_retro_ms(N_elec,:) - prop_spikes_retro_ms(1,:) );    
    
    % 2) linear regression
    x_mm = 1.0e-3 * elec_distance_um * [0 : ( N_elec - 1 ) ]'; %#ok<NBRAK>    
    fit_opt = fitoptions( 'Method', 'NonlinearLeastSquares', 'Algorithm', 'Trust-Region', 'StartPoint', [0.0, 1.0] );
    fit_fun  = fittype('x0 + v*t', 'independent', 't', 'options', fit_opt);
    
    % 2.1) antero
    N_spks        = size( prop_spikes_antero_ms, 2 );
    vel_m_per_sec = zeros( 1, N_spks );
    R2            = zeros( 1, N_spks );
    for spk = 1 : N_spks
        t_ms = prop_spikes_antero_ms(:,spk);
        [fit_obj, gof]  = fit( t_ms, x_mm, fit_fun );
        vel_m_per_sec(spk) = fit_obj.v;
        R2(spk) = gof.adjrsquare;
    end
    prop_vel_antero.regression_m_per_sec = vel_m_per_sec;
    prop_vel_antero.regression_R2        = R2;

    % 2.2) retro
    N_spks        = size( prop_spikes_retro_ms, 2 );
    vel_m_per_sec = zeros( 1, N_spks );
    R2            = zeros( 1, N_spks );
    for spk = 1 : N_spks
        t_ms = prop_spikes_retro_ms(:,spk);
        [fit_obj, gof]  = fit( t_ms, x_mm, fit_fun );
        vel_m_per_sec(spk) = fit_obj.v;
        R2(spk) = gof.adjrsquare;
    end
    prop_vel_retro.regression_m_per_sec = vel_m_per_sec;
    prop_vel_retro.regression_R2        = R2;

    
%  % VALIDATE
%     electrodes_mask = obj.Mask_From_Electrode_Labels( my_electrode_labels );
%     obj.Plot_Raster( electrodes_mask );
%     hold on
%     
%     %antero
%     for k = 1 : N_elec        
%         t = 1.0e-3 * selected_spike_times_ANT_ms(k,:);
%         y = k * ones( size(t) );
%         plot( t, y, 'go', 'MarkerSize',4 );        
%     end
%     
%     %retro
%     for k = 1 : N_elec        
%         t = 1.0e-3 * selected_spike_times_RET_ms(k,:);
%         y = k * ones( size(t) );
%         plot( t, y, 'ro', 'MarkerSize',4 );        
%     end
%        
%     antero
%     for spk = 1 : size( prop_spikes_antero_ms, 2 )
%         t = 1.0e-3 * prop_spikes_antero_ms(:,spk);
%         y = 1:N_elec;
%         plot( t, y, 'go:', 'MarkerSize',4 );        
%     end
%     
%     retro
%     for spk = 1 : size( prop_spikes_retro_ms, 2 )
%         t = 1.0e-3 * prop_spikes_retro_ms(:,spk);
%         y = 1:N_elec;
%         plot( t, y, 'ro:', 'MarkerSize',4 );        
%     end
% 
% 
%     hold off
%     title('Antero: green  |  Retro : red')
    
end