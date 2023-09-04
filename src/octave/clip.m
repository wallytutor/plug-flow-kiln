function val = clip(f_val, f_min, f_max)
    % Clip value of `f_val` to range [`f_min`, `f_max`].
    val = min(max(f_val, f_min), f_max);
end
