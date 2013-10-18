function output_struct = worst_simulink_backward(model_name, t, x,...
    nu, has_U, has_u, has_v, has_del, input_struct)


    %nu, xi, u, v, params, U, Ut, has_U, has_u, has_v, has_del)

% Get some useful variables
state_dim = size(x,2);
N = length(t);
output_dim = size(nu,2);

if has_U
    U = input_struct.U;
    Ut = input_struct.Ut;
else
    U = [];
    Ut = [];
end
if has_u
    u = input_struct.u;
    total_udim = size(u,2);
else
    u = [];
    total_udim = 0;
end
if has_v
    v = input_struct.v;
    total_vdim = size(v,2);
    xi = input_struct.xi;
else
    v = [];
    total_vdim = 0;
    xi = [];
end
if has_del
    params = input_struct.params;
    num_params = length(params);
else
    params = [];
    num_params = 0;
end




% Integrate over time - this automatically handles the fact that the system
% is a final value problem.
[time, soln] = ode45(@back_deriv, [t(N), t(1)], zeros(state_dim+num_params,1));

% Interpolate signals
soln = vector_interpolate(soln, time, t);


% Compute the output
output = zeros(N, total_udim+total_vdim);
for i = 1:N
    Ui = vector_interpolate(U, Ut, t(i))';
    if has_u, now_u = u(i,:)'; else now_u = []; end;
    if has_v, now_v = v(i,:)'; else now_v = []; end;
    [~,B,~,D] = linmod(model_name, x(i,:)', [Ui; now_u; now_v; params]);
    
    if has_u
        dfdu = B(:,1:total_udim);
        dgdu = D(1:output_dim, 1:total_udim);
        dhdu = D(output_dim+1:end, 1:total_udim);
    else
        dfdu = [];
        dgdu = [];
        dhdu = [];
    end
    if has_v
        dfdv = B(:,total_udim+1:total_udim+total_vdim);
        dgdv = D(1:output_dim, total_udim+1:total_udim+total_vdim);
        dhdv = D(output_dim+1:end, total_udim+1:total_udim+total_vdim);
    else
        dfdv = [];
        dgdv = [];
        dhdv = [];
    end
    

    Cdyn = [dfdu', zeros(total_udim, num_params); ...
            dfdv', zeros(total_vdim, num_params)];
    Ddyn = [dgdu', dhdu'; dgdv', dhdv];
    
    if has_u, now_nu = nu(i,:)'; else now_nu = []; end;
    if has_v, now_xi = xi(i,:)'; else now_xi = []; end;
    
    output(i,:) = Cdyn*soln(i,:)' + Ddyn*[now_nu; now_xi];
end
if has_u
    output_struct.gamma = -output(:, 1:total_udim);
end
if has_v
    output_struct.kappa = output(:, total_udim+1:end);
end
if has_del
    lambdad = soln(:, state_dim+1:end);
    output_struct.lambdad0 = lambdad(1,:);
end


% Nested derivative function. Uses linmod to take jacobians
function deriv = back_deriv(time, lambda)

    cur_x = vector_interpolate(x,t,time)';
    cur_u = vector_interpolate(u,t,time)';
    cur_v = vector_interpolate(v,t,time)';
    cur_nu = vector_interpolate(nu,t,time)';
    cur_xi = vector_interpolate(xi,t,time)';
    cur_U = vector_interpolate(U,Ut,time)';
    
    [Ai,Bi,Ci,Di] = linmod(model_name, cur_x, [cur_U; cur_u; cur_v; params]);
    
    dfdx = Ai;
    dgdx = Ci(1:output_dim, :);
    dhdx = Ci(output_dim+1:end, :);
    
    if has_del
        dfdd = Bi(:, total_udim+total_vdim+1:end);
        dgdd = Di(1:output_dim, total_udim+total_vdim+1:end);
        dhdd = Di(output_dim+1:end, total_udim+total_vdim+1:end);
    else
        dfdd = [];
        dgdd = [];
        dhdd = [];
    end

    Adyn = [dfdx', zeros(state_dim, num_params); ...
            dfdd', zeros(num_params, num_params)];
    Bdyn = [dgdx', dhdx'; dgdd', dhdd'];
    
    deriv = Adyn*lambda + Bdyn*[cur_nu; cur_xi];
    deriv = -deriv;
end


end