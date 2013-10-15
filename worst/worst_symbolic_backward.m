function [gamma, kappa, lambdad0] = ...
    worst_symbolic_backward(nu, xi, x, U, u, v, params, time_axis, nom_t, ...
        ...
        onestep_fdx, onestep_fdu, onestep_fdv, onestep_fddelta, onestep_hdx,...
        onestep_hdu, onestep_hdv, onestep_hddelta, onestep_gdx, onestep_gdu,...
        onestep_gdv, onestep_gddelta)

num_parameters = size(params, 1);

costate_param_dimension = num_parameters;
n = size(x, 2);
%costate_dimension = size(costate_signal, 1);
costate_dimension = n;
total_u_dim = size(u, 1);
total_U_dim = size(U, 1);
total_v_dim = size(v, 1);
time_bracket = [time_axis(1) time_axis(end)];
back_time_bracket = [time_axis(end) time_axis(1)];
N = length(time_axis);



% Compute values of dxdu, dxdv, and dxddelta, use ode45 and do this
% numerically and not symbolically
eval_x_option = 'u';
[time_dxdu, dxdu] = ode45(@eval_x_deriv, time_bracket, zeros(n, total_u_dim));
eval_x_option = 'v';
[time_dxdv, dxdv] = ode45(@eval_x_deriv, time_bracket, zeros(n, total_v_dim));
eval_x_option = 'delta';
[time_dxdd, dxdd] = ode45(@eval_x_deriv, time_bracket,zeros(n, num_parameters));







% Call the ode solver. Integrate backwards in time. 
final_condition = zeros(costate_dimension + costate_param_dimension, 1);
[solution_time, bvp_soln] = ode45(@derivative_function, back_time_bracket,...
    final_condition);
lambda = bvp_soln(:,1:n);
lambdad0 = bvp_soln(:,n+1:end);




% Interpolate the solution to the usual time axis
lambda = vector_interpolate(lambda, solution_time, time_axis);
lambdad0 = vector_interpolate(lambdad0, solution_time, time_axis);




% Compute the output and unmodeled output of the adjoint system
gamma = zeros(N, total_u_dim);
kappa = zeros(N, total_v_dim);
for j = 1:N
    current_x = x(j,:)';
    current_t = time_axis(j);
    current_u = u(j,:)';
    current_v = v(j,:)';
    current_U = vector_interpolate(U, nom_t, current_t);
    current_nu = nu(j,:)';
    current_xi = xi(j,:)';
    
    % Compute matrices
    current_dxdu = vector_interpolate(dxdu, time_dxdu, current_t);
    current_dxdv = vector_interpolate(dxdv, time_dxdv, current_t);

    symbolic_input = create_symbolic_input(current_x, current_t, current_u, ...
        current_U, current_v, params);

    dfdx = onestep_fdx(symbolic_input{:});
    dhdx = onestep_hdx(symbolic_input{:});
    dgdx = onestep_gdx(symbolic_input{:});
    dfdu = onestep_fdu(symbolic_input{:}) + dfdx*current_dxdu;
    dfdv = onestep_fdv(symbolic_input{:}) + dfdx*current_dxdv;
    dhdu = onestep_hdu(symbolic_input{:}) + dhdx*current_dxdu;
    dhdv = onestep_hdv(symbolic_input{:}) + dhdx*current_dxdv;    
    dgdu = onestep_gdu(symbolic_input{:}) + dgdx*current_dxdu;
    dgdv = onestep_gdv(symbolic_input{:}) + dgdx*current_dxdv;
    
    
    gamma(:,j) = dfdu' * lambda(:,j) + dgdu'*current_nu' + dhdu'*current_xi';
    kappa(:,j) = dfdv' * lambda(:,j) + dgdv'*current_nu' + dhdv'*current_xi';
    
    gamma = gamma * -1;
end



%------------------------------------------------------------------------------%
%                  Nested functions for the procedure above                    %
%------------------------------------------------------------------------------%


% Nested Costate derivative function
function costate_derivative = derivative_function(time, costate)  
    
    % Compue current values of dfdx and dxdu, dxdd
    cur_state = vector_interpolate(x, time_axis, time);
    cur_disturbance = vector_interpolate(u, time_axis, time);
    cur_unmodeled_input = vector_interpolate(v, ...
        time_axis, time);
    cur_nominal_input = vector_interpolate(U,time_axis,time);
    cur_dxdd = vector_interpolate(dxdd, time_dxdd, time);
    
    symb_input = create_symbolic_input(cur_state, time, cur_disturbance, ...
        cur_nominal_input, cur_unmodeled_input, params);
    
    fdx = onestep_fdx(symb_input{:});
    fdxd = onestep_fddelta(symb_input{:}) + fdx*cur_dxdd;
    hdx = onestep_hdx(symb_input{:});
    hdxd = onestep_hddelta(symb_input{:}) + hdx*cur_dxdd;
    gdx = onestep_gdx(symb_input{:});
    gdxd = onestep_gddelta(symb_input{:}) + gdx*cur_dxdd;
    
    
    % Compute the costate derivative
    cur_costate_input = vector_interpolate(nu,time_axis,time);
    cur_costate_unmodeled_input = ...
        vector_interpolate(xi, time_axis, time);
    
    costate_derivative = zeros(costate_dimension + costate_param_dimension, 1);
    
    costate_derivative(1:n) = fdx'*costate(1:n) +...
        gdx' * cur_costate_input + hdx' * cur_costate_unmodeled_input;
    costate_derivative(n+1:end) = fdxd' * ...
        costate(1:n) + ...
        gdxd' * cur_costate_input + hdxd' * cur_costate_unmodeled_input;
end


% Find time derivatives of dxdu, dxdv, and dxddelta
function deriv = eval_x_deriv(time, cur_value)
    
    cur_nominal_input = vector_interpolate(U,time_axis,time);
    cur_disturbance_input=vector_interpolate(u,time_axis,time);
    cur_unmodeled_input = vector_interpolate(v, ...
                                             time_axis, time);
                                         
    symb_input = create_symbolic_input(cur_value, time, ...
        cur_disturbance_input, cur_nominal_input, cur_unmodeled_input, ...
        params);
    
    dfdx_onestep = onestep_fdx(symb_input{:});
    if strcmp(eval_x_option, 'u')
    	offset = onestep_fdu(symb_input{:});
    elseif strcmp(eval_x_option, 'v')
    	offset = onestep_fdv(symb_input{:});        
    elseif strcmp(eval_x_option, 'delta')
    	offset = onestep_fddelta(symb_input{:});
    end
    
    deriv = offset + dfdx_onestep*cur_value;
end


function symb = create_symbolic_input(current_state, current_time, ...
        current_disturbance, current_nominal_input, ...
        current_uncertain_feedback, parameter_values)
    
    symb = [reshape(current_state, [1, n]), ...
            current_time, ...
            reshape(current_disturbance, [1, total_u_dim]), ...
            reshape(current_nominal_input, [1, total_U_dim]), ...
            reshape(current_uncertain_feedback, ...
                    [1,total_v_dim]),...
            reshape(parameter_values, [1, num_parameters])];
    symb = num2cell(symb);
end

end