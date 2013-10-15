function state_signal = worst_symbolic_foward(...
    state_derivative, initial_condition, parameters, disturbance_signal, ...
    unmodeled_feedback_input, nominal_input, time_axis, nominal_time_axis)

% Use ode45 to solve the IVP
times = [time_axis(1) time_axis(end)];
[solution_time, state_signal] = ode45(@df, times, initial_condition);
state_signal = vector_interpolate(state_signal, solution_time, time_axis);

% Nested derivative function
function derivative = df(time, state)

    current_u = vector_interpolate(disturbance_signal,time_axis,time);
    current_U = vector_interpolate(nominal_input, nominal_time_axis, time);
    current_v = vector_interpolate(unmodeled_feedback_input, time_axis, time);

    % Standard non-augmented state derivative
    derivative = state_derivative(state, time, current_u, current_U, ...
        current_v, parameters);

end

end