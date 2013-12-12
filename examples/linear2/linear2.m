% This example computes the worst possible disturbance of a linear system
% with the following dynamics:
%
%     x1'(t) = -3*x1(t) + a*x2(t) + u(t)
%     x2'(t) = 2*x1(t) - 5*x2(t) + 2*u(t)
%     y(t)   = x1(t)
%
% where (x1, x2)' is the state, u(t) is a one-dimensional disturbance of norm 
% 1, and y(t) is the output.



clear; clc; close all;


model_name = 'linear2_model';

disturbance_specs = [1 1];
params = [-1 0 1];

ti_val = 0; 
tf_val = 10;

output_dim = 1;

error_tol = .01;


output_struct = ...
    worst(model_name, output_dim, 'ti', ti_val, 'tf', tf_val, ...
          'disturbance_specs', disturbance_specs, 'error_tol', error_tol, ...
          'params', params, 'num_iter', 1);
      
% Plot the distribution of worst values
figure
hist(output_struct.costs)