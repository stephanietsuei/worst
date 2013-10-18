% This example computes the worst possible disturbance of a linear system
% with the following dynamics:
%
%     x1'(t) = -3*x1(t) + a*x2(t) + u(t) + 2*v(t)
%     x2'(t) = 2*x1(t) - 5*x2(t) + 2*u(t) + v(t)
%     y(t)   = x1(t) + v(t)
%     z(t)   = -x2(t) + v(t)
%
% where (x1, x2)' is the state, u(t) is a one-dimensional disturbance of norm 
% 1, y(t) is the output, and v(t) and z(t) form a feedback loop that models
% uncertain dynamics; v(t) and z(t) have the same finite-time 2-norm. 



clear; clc; close all;


model_name = 'linear1_model';

disturbance_specs = [1 1];
unmodeled_io = [1 1];
params = [-1 0 1];

ti_val = 0; 
tf_val = 10;

output_dim = 1;

max_iterations = 10;
error_tol = .01;


output_struct = ...
    worst('simulink', model_name, output_dim, 'ti', ti_val, 'tf', tf_val, ...
          'disturbance_specs', disturbance_specs, 'error_tol', error_tol, ...
          'unmodeled_io', unmodeled_io, 'params', params);
      
display(['Total cost is: ' num2str(output_struct.cost)]);

figure
plot(output_struct.time_axis, output_struct.d)
title(['Worst possible disturbance'])
xlabel('Time')
ylabel('Disturbance')