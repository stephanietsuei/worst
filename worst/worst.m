function output_struct = worst(varargin)


allowed_modes = {'symbolic', 'simulink'};
assert(ismember(varargin{1}, allowed_modes))



% Parser for simulink input
simulink_parse = inputParser;
addRequired(simulink_parse, 'system_name', @(x) isa(x,'char'));
addRequired(simulink_parse, 'output_dim', @(x) isa(x,'double'));
addParamValue(simulink_parse, 'tf', 10, @(x) isa(x,'double'));
addParamValue(simulink_parse, 'ti', 0, @(x) isa(x, 'double'));
addParamValue(simulink_parse, 'nominal_input', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'nominal_time', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'disturbance_specs', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'unmodeled_io', [], @(x) isa(x,'double'));
addParamValue(simulink_parse, 'params', [], @(x) isa(x, 'double'));
addParamValue(simulink_parse, 'max_iter', 40, @(x) isa(x,'double'));
addParamValue(simulink_parse, 'error_tol', .01, @(x) isa(x,'double'));


% Parser for symbolic input
symbolic_parse = inputParser;
% TODO


% Parse input for each type of input
if strcmp(varargin{1}, 'symbolic')
    parse(symbolic_parse, varargin{2:end});
elseif strcmp(varargin{1}, 'simulink')
    parse(simulink_parse, varargin{2:end});
    disturbance_specs = simulink_parse.Results.disturbance_specs;
    error_tol = simulink_parse.Results.error_tol;
    max_iter = simulink_parse.Results.max_iter;
    nominal_input = simulink_parse.Results.nominal_input;
    nominal_time = simulink_parse.Results.nominal_time;
    output_dim = simulink_parse.Results.output_dim;
    params = simulink_parse.Results.params;
    system_name = simulink_parse.Results.system_name;
    tf = simulink_parse.Results.tf;
    ti = simulink_parse.Results.ti;
    unmodeled_io = simulink_parse.Results.unmodeled_io;
end


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


% Call worst
if strcmp(varargin{1}, 'simulink')
    output_struct = worst_simulink(system_name, disturbance_specs, ...
        unmodeled_io, params, nominal_input, nominal_time, ti, tf, max_iter, ...
        output_dim, error_tol);
elseif strcmp(varargin{1}, 'symbolic')
    output_struct = worst_symbolic();
end



end