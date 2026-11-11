function f = check(name, got, expect, rtol)
err = abs(got - expect) / max(abs(expect), eps);
f = err > rtol;
if f, tag = 'FAIL'; else, tag = 'pass'; end
fprintf('%s  %-40s got %.6g  expected %.6g  (rel err %.1e)\n', tag, name, got, expect, err);
end
