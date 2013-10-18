% Computes the inner product of two multidimensional time signals. 
% The two signals should be formatted as such:
%   u = [u1; u2; u3; ... un]
% where the rows ui are
%   ui = [ui(t1) ui(t2) ui(t3) ... ui(tn)]
% The time axis should be a row vector

function value = multidim_ip(signal1, signal2)

    % Sanity check - make sure that signals are compatible
    assert(isequal(size(signal1), size(signal2)));
    
    block = signal1' * signal2;
    value = sum(diag(block));

end