function value = multidim_norm(signal, time_axis)

time_range = time_axis(end) - time_axis(1);
value = sqrt(sum(trapz(time_axis, signal.^2)) / time_range);

end