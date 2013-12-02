function out = worst_simulink(model_name, disturbance_specs, unmodeled_io, ...
    params, nominal_input, nominal_t, ti, tf, max_iterations, output_dim, ...
    error_tol, averaging, nominal_output)

% Computes the worst case error norm of a system given a proper simulink
% model. For now, this file assumes the following


% Input order of simulink model
% 1. nominal input (in one block)
% 2. disturbances are next input (in one block)
% 3. unmodeled inputs are next input (in one block)
% 4. uncertain parameters (in one block)

% Output order of simulink model
% 1. error signal
% 2. unmodeled feedback model

% Input format of "disturbances". Each row represents a separate signal
% [dimension, norm]

% Input format of unmodeled io
% [input dimension, output dimension]

% input format of parameters
% [lower bound, nominal value, upper bound]

% ti and tf are scalars where tf > ti

% assumes initial state is either specified in the model or is located in
% the base workspace.





%------------------------------------------------------------------------------%
% Determine what inputs we have and get some useful variables
%------------------------------------------------------------------------------%
has_U = 1;
has_u = 1;
has_v = 1;
has_del = 1;
using_out = 0;

if isempty(disturbance_specs), has_u = 0; end
if isempty(unmodeled_io), has_v = 0; end
if isempty(params), has_del = 0; end
if isequal([has_u,has_v,has_del], [0 0 0])
    display('Why are you even using this program?')
    exit;
end

if isempty(nominal_input)
    has_U = 0;
end
if ~isempty(nominal_output)
    using_out = 1;
end

simulink_input_str1 = '';
simulink_input_str2 = '';
if has_U
    simulink_input_str1 = 'nomnom_U';
    simulink_input_str2 = 'nomnom_U';
end
if has_u
    simulink_input_str1 = [simulink_input_str1, ',', 'nomnom_u1'];
    simulink_input_str2 = [simulink_input_str2, ',', 'nom_u'];
end
if has_v
    simulink_input_str1 = [simulink_input_str1, ',', 'nomnom_v1'];
    simulink_input_str2 = [simulink_input_str2, ',', 'nom_v'];
end
if has_del
    simulink_input_str1 = [simulink_input_str1, ',', 'nomnom_del1'];
    simulink_input_str2 = [simulink_input_str2, ',', 'nom_delta'];
end


% Get useful variables
if has_u, 
    total_disturbance_dim = sum(disturbance_specs(:,1));
else
    total_disturbance_dim = 0;
end

num_unmodeled_in = size(unmodeled_io,1);
if has_v
    total_v_dim = sum(unmodeled_io(:,1));
    total_z_dim = sum(unmodeled_io(:,2));
else
    total_v_dim = 0;
    total_z_dim = 0;
end

num_params = size(params,1);
if has_del
    pars = params(:,2)';
else
    pars = [];
end



%------------------------------------------------------------------------------%
% First, simulate the system with zero distubances to get a time axis and a
% nominal output.
%------------------------------------------------------------------------------%

% Set Simulink simulation parameters
sim_parameters1.StartTime = 'nomnom_ti';
sim_parameters1.StopTime = 'nomnom_tf';
sim_parameters1.LoadExternalInput = 'on';
sim_parameters1.ExternalInput = simulink_input_str1;
assignin('base', 'nomnom_ti', ti);
assignin('base', 'nomnom_tf', tf);
    
% Make empty inputs, and then put the inputs in the base workspace where
% Simulink can find them
if has_u
    u1.time = [ti; tf];
    u1.signals.values = zeros(2, total_disturbance_dim);
    u1.signals.dimensions = total_disturbance_dim;
    assignin('base', 'nomnom_u1', u1);
end
if has_v
    v1.time = [ti; tf];
    v1.signals.values = zeros(2, total_v_dim);
    v1.signals.dimensions = total_v_dim;
    assignin('base', 'nomnom_v1', v1);
end
if has_del
    del1.time = [ti; tf];
    del1.signals.values = [pars; pars];
    del1.signals.dimensions = num_params;
    assignin('base', 'nomnom_del1', del1);
end
if has_U
    U.time = nominal_t;
    U.signals.values = nominal_input;
    U.signals.dimensions = size(nominal_input,2);
    assignin('base', 'nomnom_U', U);
end
    
% Simulate the system if we weren't given a nominal output
if ~using_out
    
    % Set simulation parameters and simulate the system if we weren't given a
    % nominal output
    simout = sim(model_name, sim_parameters1);
    time_axis = simout.get('tout');
    nominal_output = simout.get('yout');
    nominal_output = nominal_output(:,1:output_dim);
else
    time_axis = nominal_t;
end

N = length(time_axis);


%------------------------------------------------------------------------------%
% Initialize worst case signals
%------------------------------------------------------------------------------%

% Initialize disturbance and unmodeled input signals, normalize
% disturbance, and then initialize parameters input
if has_u
    u.time = time_axis;
    u.signals.dimensions = total_disturbance_dim;
    u.signals.values = rand(N, sum(disturbance_specs(:,1)));
    u.signals.values = normalize_signal(u.signals.values, disturbance_specs, ...
        time_axis);
end
if has_v
    v.time = time_axis;
    v.signals.values = zeros(N, sum(unmodeled_io(:,1)));
    v.signals.dimensions = total_disturbance_dim;
end
if has_del
    delta.time = time_axis;
    delta.signals.values = repmat(params(:,2)', N, 1);
    delta.signals.dimensions = num_params;
end


% Initialize costate for unmodeled inputs
if has_v, lambdaD = ones(1, num_unmodeled_in); end


%------------------------------------------------------------------------------%
% The iteration
%------------------------------------------------------------------------------%

% Start the iteration
converged = 0;
iterations = 0;
output = Inf*ones(N, output_dim);

options_struct = sim_parameters1;
options_struct.ExternalInput = simulink_input_str2;
options_struct.SaveState = 'on';


while ((~converged) && (iterations <= max_iterations))
    
    display(iterations);
    
    last_output = output;
    
    
    % Integrate the system forward in time
    if has_u, assignin('base','nom_u',u); end
    if has_v, assignin('base','nom_v',v); end
    if has_del, assignin('base','nom_delta',delta); end
    simout = sim(model_name, options_struct);
    t = simout.get('tout');
    yout = simout.get('yout');
    output = vector_interpolate(yout(:,1:output_dim), t, time_axis);
    output = output - nominal_output;
    if has_v
        z = vector_interpolate(yout(:,output_dim+1:end), t, time_axis);
    end
    x = vector_interpolate(simout.get('xout'), t, time_axis);
    
    
    % Align the system and get inputs for third part
    nu = output;
    if has_U
        back_struct.U = U.signals.values;
        back_struct.Ut = U.time;
    end
    if has_del
        back_struct.params = delta.signals.values(1,:)';
    end
    if has_u
        back_struct.u = u.signals.values;
    end
    if has_v
        start_ind = 1;
        xi = zeros(N, total_z_dim);
        for i = 1:num_unmodeled_in
            end_ind = start_ind + unmodeled_io(i,2) - 1;
            xi(:, start_ind:end_ind) = lambdaD(i) * z(:, start_ind:end_ind);
            start_ind = end_ind + 1;
        end
        back_struct.xi = xi;
        back_struct.v = v.signals.values;
    end
    
    
    % Integrate backward in time
    back_output = worst_simulink_backward(model_name, time_axis, x, nu, ...
        has_U, has_u, has_v, has_del, back_struct);
    if has_u, gamma = back_output.gamma; end
    if has_v, kappa = back_output.kappa; end
    if has_del, lambdad0 = back_output.lambdad0; end
    
    
    % Realign the system
    if has_u
        realign_input.gamma = gamma;
        realign_input.u_specs = disturbance_specs;
        realign_input.u = u.signals.values;
    end
    if has_v
        realign_input.kappa = kappa;
        realign_input.io_specs = unmodeled_io;
        realign_input.v = v.signals.values;
        realign_input.z = z;
    end
    if has_del
        realign_input.lambdad0 = lambdad0;
        realign_input.delta = delta.signals.values(1,:)';
        realign_input.delta = [params(:,1) realign_input.delta params(:,3)];
    end
    realign_output = worst_realign(time_axis, has_u, has_v, has_del, ...
        realign_input); 

    
    % Store new values for next iteration
    if has_u
        old_u = u.signals.values;
        u.signals.values = realign_output.u; 
    end
    if has_v
        old_v = v.signals.values;
        v.signals.values = realign_output.v;
        lambdaD = realign_output.lambdaD;
    end
    if has_del, delta.signals.values = repmat(realign_output.params', N, 1); end
    
    
    % If we're averaging values together for stability, then repeat the
    % realign step with the new values of u, v, and z
    if averaging
        if has_u, realign_input.u = .5*(u.signals.values + old_u); end
        if has_v, realign_input.v = .5*(v.signals.values + old_v); end
        if (has_u || has_v)
            realign_output = worst_realign(time_axis, has_u, has_v, 0, ...
                realign_input);
        end
        if has_u, u.signals.values = realign_output.u; end
        if has_v,
            v.signals.values = realign_output.v;
            lambdaD = realign_output.lambdaD;
        end
    end
    
    
    % Compute whether or not we're done
    error = max(sum((last_output-output).^2)./sum(output.^2));
    if (error < error_tol)
        converged = 1;
    end
    iterations = iterations + 1;
    display(error)
end





%------------------------------------------------------------------------------%
% Compute final cost and get final values
%------------------------------------------------------------------------------%

if has_u, out.d = u.signals.values; end;
if has_v, out.v = v.signals.values; end;
if has_del, out.parm = delta.signals.values(1,:)'; end;
out.time_axis = time_axis;
out.cost = sum(trapz(time_axis, (last_output-output).^2));
out.converged = converged;

end