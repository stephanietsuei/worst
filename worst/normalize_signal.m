function normalized_signals = normalize_signal(signal, signal_specs, time_axis)


subsignal_dims = signal_specs(:,1);
subsignal_norms = signal_specs(:,2);
num_signals = length(subsignal_dims);



% get start indices in signal columns
start_indices = zeros(num_signals,1);
start_indices(1) = 1;
for i = 2:1:num_signals  % this won't run if there aren't enough dimensions
    start_indices(i) = sum(subsignal_dims(1:i-1)) + 1;
end


% normalize signals
normalized_signals = zeros(size(signal));
for i = 1:num_signals
    start_ind = start_indices(i);
    end_ind = start_ind + subsignal_dims(i) - 1;
    
    current_signal = signal(:,start_ind:end_ind);
    current_norm = multidim_norm(current_signal, time_axis);
    
    normalized_signals(:,start_ind:end_ind) = ...
        current_signal * subsignal_norms(i) / current_norm;
end


end