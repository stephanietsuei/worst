function output = worst(varargin)


% Parser for simulink input
simulink_parse = inputParser;
addRequired(simulink_parse, 'system_name', @(x) isa(x,'char'));
addRequired(simulink_parse, 'output_dim', @(x) isa(x,'double'));
addParamValue(simulink_parse, 'tf', 10, @(x) isa(x,'double'));
addParamValue(simulink_parse, 'ti', 0, @(x) isa(x, 'double'));
addParamValue(simulink_parse, 'nominal_input', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'nominal_output', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'nominal_time', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'disturbance_specs', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'unmodeled_io', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'params', [], @(x) isa(x, 'double'));
addParamValue(simulink_parse, 'max_iter', 40, @(x) isa(x,'double'));
addParamValue(simulink_parse, 'error_tol', .01, @(x) isa(x,'double'));
addParamValue(simulink_parse, 'num_iter', 10, @(x) ((mod(x,1)==0) && (x > 0)))
addParamValue(simulink_parse, 'averaging', 1, @(x) ismember(x, [0 1]))
addParamValue(simulink_parse, 'plot_cost', 0, @(x) ismember(x, [0 1]))
addParamValue(simulink_parse, 'plot_d', 0, @(x) ismember(x, [0 1]))
addParamValue(simulink_parse, 'plot_v', 0, @(x) ismember(x, [0 1]))
addParamValue(simulink_parse, 'plot_parm', 0, @(x) ismember(x, [0 1]))
addParamValue(simulink_parse, 'plot_error', 0, @(x) ismember(x, [0 1]))


% Parse input for each type of input
parse(simulink_parse, varargin{:});
disturbance_specs = simulink_parse.Results.disturbance_specs;
error_tol = simulink_parse.Results.error_tol;
max_iter = simulink_parse.Results.max_iter;
nominal_input = simulink_parse.Results.nominal_input;
nominal_output = simulink_parse.Results.nominal_output;
nominal_time = simulink_parse.Results.nominal_time;
output_dim = simulink_parse.Results.output_dim;
params = simulink_parse.Results.params;
system_name = simulink_parse.Results.system_name;
tf = simulink_parse.Results.tf;
ti = simulink_parse.Results.ti;
unmodeled_io = simulink_parse.Results.unmodeled_io;
num_iter = simulink_parse.Results.num_iter;
averaging = simulink_parse.Results.averaging;
plot_cost = simulink_parse.Results.plot_cost;
plot_d = simulink_parse.Results.plot_d;
plot_v = simulink_parse.Results.plot_v;
plot_parm = simulink_parse.Results.plot_parm;
plot_error = simulink_parse.Results.plot_error;


% Some more input checking
assert(tf > ti);
assert(max_iter > 1);
assert((0 < error_tol) & (error_tol < 1));
if ~isempty(nominal_time)
    assert(ti == nominal_time(1))
    assert(tf == nominal_time(end))
end
if ~isempty(disturbance_specs)
    assert(size(disturbance_specs, 2) == 2)
    assert(isempty(find(disturbance_specs<=0, 1)))
end
if ~isempty(unmodeled_io)
    assert(size(unmodeled_io, 2) == 2)
    assert(isempty(find(unmodeled_io<=0, 1)))
end
if ~isempty(params)
    assert(size(params, 2) == 3);
    num_params = size(params, 1);
    for i = 1:num_params
        assert((params(i,1) <= params(i,2)) & (params(i,2) <= params(i,3)))
        assert(params(i,1) < params(i,3))
    end
end


% Call worst num_iter times and store its output
output.costs = zeros(1,num_iter);
output.converged = zeros(1,num_iter);
output.time_axis = cell(1,num_iter);
if ~isempty(params), output.parm = cell(1,num_iter); end;
if ~isempty(disturbance_specs), output.d = cell(1,num_iter); end;
if ~isempty(unmodeled_io), output.v = cell(1,num_iter); end;
for i = 1:num_iter
    output_struct = worst_simulink(system_name, disturbance_specs, ...
        unmodeled_io, params, nominal_input, nominal_time, ti, tf, max_iter, ...
        output_dim, error_tol, averaging, nominal_output, plot_cost, plot_d, ...
        plot_v, plot_parm, plot_error);
    output.costs(i) = output_struct.cost;
    output.converged(i) = output_struct.converged;
    output.time_axis{i} = output_struct.time_axis;
    if ~isempty(params), output.parm{i} = output_struct.parm; end;
    if ~isempty(disturbance_specs), output.d{i} = output_struct.d; end;
    if ~isempty(unmodeled_io), output.v{i} = output_struct.v; end;
end


end