function output_struct = worst_realign(time_axis, has_u, has_v, has_del, ...
    input_struct)

N = length(time_axis);

if has_u
    disturbance_specs = input_struct.u_specs;
    disturbance_dims = disturbance_specs(:,1);
    disturbance_norms = disturbance_specs(:,2);
    num_disturbances = length(disturbance_dims);
    
    u = input_struct.u;
    gamma = input_struct.gamma;

    lambdau = zeros(num_disturbances, 1);
    new_u = zeros(N, sum(disturbance_dims));
end

if has_v
    unmodeled_io = input_struct.io_specs;
    num_unmodeled_loops = size(unmodeled_io,1);
    
    v = input_struct.v;
    z = input_struct.z;
    kappa = input_struct.kappa;
    
    lambdaD = zeros(num_unmodeled_loops, 1);
    new_v = zeros(N, sum(unmodeled_io(:,1)));
end

if has_del
    params = input_struct.delta;
    lambdad0 = input_struct.lambdad0;
end



% Compute norms of z signals
if has_v
    z_norms = zeros(num_unmodeled_loops, 1);
    start_ind = 1;
    for i = 1:num_unmodeled_loops
        end_ind = start_ind + unmodeled_io(i,2) - 1;
        z_norms(i) = multidim_norm(z(:,start_ind:end_ind), time_axis);
        start_ind = end_ind + 1;
    end
end


% Compute lambdau and u
if has_u
    start_ind = 1;
    for i = 1:num_disturbances
        end_ind = start_ind + disturbance_dims(i) - 1;

        cur_gam = gamma(:, start_ind:end_ind);
        cur_u = u(:, start_ind:end_ind);
        cur_gam_norm = multidim_norm(cur_gam, time_axis);
        signs = sign(multidim_ip(cur_gam, cur_u));
        if signs == 0
            signs = 1;
        end

        lambdau(i) = signs * cur_gam_norm / disturbance_norms(i);
        new_u(:, start_ind:end_ind) = cur_gam / lambdau(i);

        start_ind = end_ind + 1;
    end
    output_struct.u = new_u;
end


% Compute lambdaD and v
if has_v
    start_ind = 1;
    for i = 1:num_unmodeled_loops
        end_ind = start_ind + unmodeled_io(i,1) - 1;

        cur_kap = kappa(:, start_ind:end_ind);
        cur_v = v(:, start_ind:end_ind);
        cur_kap_norm = multidim_norm(cur_kap, time_axis);
        signs = sign(multidim_ip(cur_kap, cur_v));
        if signs == 0
            signs = 1;
        end

        lambdaD(i) = signs * cur_kap_norm / z_norms(i);
        new_v(:, start_ind:end_ind) = cur_kap / lambdaD(i);

        start_ind = end_ind + 1;
    end
    output_struct.v = new_v;
    output_struct.lambdaD = lambdaD;
end

% Compute new parameter values
if has_del
    param_list = params(:,2);
    new_params = param_list + lambdad0;
    new_params = min([new_params, params(:,3)]);
    new_params = max([new_params, params(:,1)]);
    output_struct.params = new_params;
end



end