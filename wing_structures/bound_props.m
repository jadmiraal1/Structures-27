function p = bound_props(m, k)
% BOUND_PROPS  Collapse a material struct with [low high] fields to scalars.
%   k = 1 -> low bound, k = 2 -> high bound. Scalar fields pass through.
p = m;
f = fieldnames(m);
for i = 1:numel(f)
    v = m.(f{i});
    if isnumeric(v) && numel(v) == 2
        p.(f{i}) = v(k);
    end
end
end
