function [f] = flatten(M)
    % Flatten matrix row-wise columnar format.
    f = reshape(M', numel(M), 1);
end