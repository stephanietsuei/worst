function [cost, u, parameters, error_signal, ...
          v, converge] = ...
   worst_symbolic(state_derivative, output_function, uncertainty_map, ...
                  initial_condition, nominal_time_axis, nominal_input, ...
                  output_dim, disturbance_definitions, ...
                  unmodeled_io, parameter_matrix, ...
                  error_tol, max_iter)



% Useful variables
time_span = [nominal_time_axis(1), nominal_time_axis(end)];

state_dimension = length(initial_condition);

%num_disturbances = length(disturbance_definitions(:,1));
total_u_dim = sum(disturbance_definitions(:,1));
%disturbance_norms = disturbance_definitions(:,2);

num_parameters = size(parameter_matrix, 1);
parameters = parameter_matrix(:,2);

total_v_dim = sum(unmodeled_io(:,1));
total_z_dim = sum(unmodeled_io(:,2));
num_unmodeled_blocks = length(unmodeled_io(:,1));

total_U_dim = size(nominal_input, 1);



%------------------------------------------------------------------------------%
% Compute the nominal trajectory, get a time axis
%------------------------------------------------------------------------------%

% Fake signals used to generate the nominal output
short_disturb_input = zeros(total_u_dim, 1);
short_unmodeled_input = zeros(total_v_dim, 1);
[time_axis, x] = ode45(@nominal_traj_deriv, time_span, ...
    initial_condition);
N = length(time_axis);

% Inputs for generating the nominal output trajectory
u = zeros(N, total_u_dim);
v = zeros(N, total_v_dim);

% Generate the nominal output
nominal_output = eval_signal(time_axis, output_function, output_dim);

% Derivative for integrating ode45
function deriv = nominal_traj_deriv(time, state)
    input = vector_interpolate(nominal_input, nominal_time_axis, time);
    deriv = state_derivative(state, time, short_disturb_input, input, ...
                             short_unmodeled_input, parameter_matrix(:,2));
end



%------------------------------------------------------------------------------%
% Initialize Variables for the power iteration
%------------------------------------------------------------------------------%


% Change the value of u for the iteration (initial value of v stays zero)
u = randn(N, total_u_dim);
u = normalize_signal(u, disturbance_definitions, time_axis);

% Iteration control
error = Inf;
error_signal = Inf*ones(N, output_dim);
num_iterations = 0;

% Costate Variables
lambdaD = ones(1, num_unmodeled_blocks);




%------------------------------------------------------------------------------%
% Generate the symbolic derivative functions
%------------------------------------------------------------------------------%

% Create symbolic variables. These are not signals, they are all row
% vectors. (Not column vectors because MATLAB 
sym_x = sym('sym_x', [1 state_dimension]);
sym_t = sym('sym_t', [1 1]);
sym_u = sym('sym_u', [1 total_u_dim]);
sym_U = sym('sym_U', [1 total_U_dim]);
sym_v = sym('sym_v', [1 total_v_dim]);
sym_delta = sym('sym_delta', [1 num_parameters]);

symbolic_input_vector = [sym_x, sym_t, sym_u, sym_U, sym_v, sym_delta];

% Convert numerical closed-form functions into symbolic functions and take
% their derivatives
f = symfun(state_derivative(sym_x, sym_t, sym_u, sym_U, sym_v, sym_delta),...
           symbolic_input_vector);
h = symfun(uncertainty_map(sym_x, sym_t, sym_u, sym_U, sym_v, sym_delta),...
           symbolic_input_vector);
Y = symfun(output_function(sym_x, sym_t, sym_u, sym_U, sym_v, sym_delta),...
           symbolic_input_vector);

% Compute half-derivatives of f, g, and h. These returned derivatives are 
% callable functions where you can just put in the current time and inputs
% and get a return value
onestep_fdx = matlabFunction(jacobian(f, sym_x));
onestep_fdu = matlabFunction(jacobian(f, sym_u));
onestep_fdv = matlabFunction(jacobian(f, sym_v));
onestep_fddelta = matlabFunction(jacobian(f, sym_delta));
onestep_hdx = matlabFunction(jacobian(h, sym_x));
onestep_hdu = matlabFunction(jacobian(h, sym_u));
onestep_hdv = matlabFunction(jacobian(h, sym_v));
onestep_hddelta = matlabFunction(jacobian(h, sym_delta));
onestep_gdx = matlabFunction(jacobian(Y, sym_x));
onestep_gdu = matlabFunction(jacobian(Y, sym_u));
onestep_gdv = matlabFunction(jacobian(Y, sym_v));
onestep_gddelta = matlabFunction(jacobian(Y, sym_delta));






%------------------------------------------------------------------------------%
% The Iteration 
%------------------------------------------------------------------------------%
while (error >= error_tol && num_iterations <= max_iter)
    
    previous_err_signal = error_signal;
    display(num_iterations);
    
    
    % FIRST MAP: Integrate system forward in time to get input and output
    x = worst_symbolic_foward(state_derivative, initial_condition,...
        parameters, u, v, nominal_input, time_axis, nominal_time_axis);
    output_signal = eval_signal(time_axis, output_function, output_dim);
    error_signal = output_signal - nominal_output;
    z = eval_signal(time_axis, uncertainty_map, total_z_dim);
    
    
    
    % SECOND MAP: Align System
    xi = zeros(N, total_z_dim);
    start_ind = 1;
    for i = 1:num_unmodeled_blocks
        end_ind = start_ind + unmodeled_io(i,2) - 1;
        xi(:,start_ind:end_ind) = z(:,start_ind:end_ind) * lambdaD(i);
        start_ind = end_ind + 1;
    end     
    nu = error_signal;
        %{
    repeated_lambda_Delta = repeat_entries(costate_unmodeled,unmodeled_io(:,2));
    costate_input = output_signal;
    costate_unmodeled_input = vector_signal_multiplication(...
         repeated_lambda_Delta, unmodeled_feedback_output);
        %}
    
    
    
    % THIRD MAP: Integrate adjoint system backward in time
    [gamma, kappa, lambdad0] = worst_symbolic_backward(nu, xi, x, ...
        nominal_input, u, v, parameters, time_axis, ...
        ...
        onestep_fdx, onestep_fdu, onestep_fdv, onestep_fddelta, onestep_hdx,...
        onestep_hdu, onestep_hdv, onestep_hddelta, onestep_gdx, onestep_gdu,...
        onestep_gdv, onestep_gddelta);
    
    
    
    % FOURTH MAP: ALIGN SYSTEM AND SET NEW VALUES
    [lambdaD, u, v, parameters] = worst_realign(gamma, kappa, parameters, ...
        lambdad0, disturbance_definitions, unmodeled_io, u, v, z, time_axis);
    %{
    % Compute the norms of each of the unmodeled signal blocks
    sizes = zeros(num_unmodeled_blocks, 1);
    start_index = 1;
    for i = 1:num_unmodeled_blocks
        end_index = start_index + unmodeled_io(i,2) - 1;
        signal = z(start_index:end_index,:);
        signal_norms = sqrt(trapz(time_axis, (signal.^2)'));
        sizes(i) = sum(signal_norms);
        start_index = end_index + 1;
    end
    
    % new costate unmodeled dynamics input (lambda_Delta)
    costate_unmodeled = zeros(num_unmodeled_blocks, 1);
    numerator = v .* v;
    %numerator = costate_unmodeled_output .* unmodeled_feedback_input;
    numerator = trapz(time_axis, numerator')';
    start_index = 1;
    for i = 1:num_unmodeled_blocks
        end_index = start_index + unmodeled_io(i,1) - 1;
        signal_norms = numerator(start_index:end_index);
        %num = sign(sum(numerator2(start_index:end_index)));
        %if sizes(i) ~= 0
            costate_unmodeled(i) = sum(signal_norms) / sizes(i)^2;
        %else
        %    costate_unmodeled(i) = sum(signal_norms);
        %end
        %costate_unmodeled(i) = costate_unmodeled(i)*num;
        %costate_unmodeled(i) = costate_unmodeled(i)*sign(costate_unmodeled(i));
        start_index = end_index + 1;
    end
    

    % new costate disturbance state
    costate_disturbance = zeros(num_disturbances, 1);
    numerator = u .* u;
    numerator = trapz(time_axis, numerator')';
    start_index = 1;
    for i = 1:num_disturbances
        end_index = start_index + disturbance_definitions(i,1) - 1;
        signal_norms = numerator(start_index:end_index);
        costate_disturbance(i) = sum(signal_norms)/disturbance_norms(i)^2;
        start_index = end_index + 1;
    end
    
    
    % new disturbance input
    repeated_vector = repeat_entries(costate_disturbance, ...
        disturbance_definitions(:,1));
    repeated_vector = elementwise_reciprocal(repeated_vector);
    u = vector_signal_multiplication(repeated_vector, ...
        costate_output);
    u = normalize_signal(u, ...
        disturbance_norms, disturbance_definitions(:,1));
    
    % new unmodeled feedback input (v)
    repeated_vector = repeat_entries(costate_unmodeled, ...
        unmodeled_io(:,1));
    repeated_vector = elementwise_reciprocal(repeated_vector);
    v = vector_signal_multiplication(...
        repeated_vector, costate_unmodeled_output);
    v = normalize_signal(v, ...
        sizes, unmodeled_io(:,1));
    
    % new nominal parameter values
    for i = 1:length(parameters)
        temp_var = parameters + costate_parameters(i,1);
        if temp_var < parameter_matrix(i,1)
            parameters(i) = parameter_matrix(i,1);
        elseif temp_var > parameter_matrix(i,3)
            parameters(i) = parameter_matrix(i,3);
        else
            parameters(i) = temp_var;
        end
    end
    %}
    
    
    % END OF ITERATION, RECOMPUTE ERRORS, SET PREVIOUS SIGNALS TO NEW
    % VALUES
    error = max(sum((error_signal - previous_err_signal).^2) ./ ...
                sum((error_signal).^2));
    display(error)
    num_iterations = num_iterations + 1;
end


% Did we converge?
if num_iterations <= max_iter
    converge = true;
else
    converge = false;
end

% Compute the cost
cost = trapz(time_axis, error_signal.^2);




%------------------------------------------------------------------------------%
% Useful nested functions 
%------------------------------------------------------------------------------%

%{
% Normalizing the disturbance
function normed = normalize_signal(signal, norms, block_dimensions)
    num_blocks = length(norms);
    normed = zeros(size(signal));
    start_ind = 1;
    for ind = 1:num_blocks
        end_ind = start_ind + block_dimensions(1) - 1;
        signal_piece = signal(start_ind:end_ind,:);
        signal_norm = sqrt(trapz(time_axis, (signal_piece.^2)'));
        signal_norm = sum(signal_norm,2);
        ratio = signal_norm / norms(ind);
        normed(start_ind:end_ind,:) = signal_piece / ratio;
        start_ind = end_ind + 1;
    end
end
%}


% Evaluating signals
function signal = eval_signal(time_axis, inst_function, signal_dimension)
    signal = zeros(length(time_axis), signal_dimension);
    for ind = 1:length(time_axis)
        signal(ind,:) = inst_function(x(:,ind), time_axis(ind), ...
            u(:,ind), nominal_input(:,ind), v(:,ind), parameters);
    end
end


%{
% Nested function for vector-signal multiplication
function new_signal = vector_signal_multiplication(vector, signal)
    dimension = size(vector, 1);
    assert(dimension == size(signal, 1))
    new_signal = zeros(size(signal));
    for dim = 1:dimension
        new_signal(dim,:) = vector(dim) * signal(dim,:);
    end
end

% Take the reciprocal of every element in a vector
function new_vector = elementwise_reciprocal(vector)
    new_vector = ones(size(vector)) ./ vector;
end


% Repeat vectors
function repeated = repeat_entries(vector, times_to_repeat)
    assert(length(vector) == length(times_to_repeat))
    assert(~length(find(times_to_repeat==0)))
    
    repeated = zeros(sum(times_to_repeat), 1);
    cur_start = 1;
    for ind = 1:length(vector)
        cur_end = cur_start + times_to_repeat(ind) - 1;
        repeated(cur_start:cur_end) = vector(ind)*ones(cur_end-cur_start+1,1);
        cur_start = cur_end + 1;
    end
end
%}

end