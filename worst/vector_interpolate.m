function new_vector = vector_interpolate(old_vector, old_time, new_time)

dimension = size(old_vector,2);
new_vector = zeros(length(new_time), dimension);

for i = 1:dimension
    new_vector(:,i) = interp1(old_time, old_vector(:,i), new_time);
end

end