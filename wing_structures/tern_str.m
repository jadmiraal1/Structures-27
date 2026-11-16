function s = tern_str(v)
% TERN_STR  Value to CSV text (numbers, logicals, strings).
if ischar(v), s = ['"' v '"'];
elseif islogical(v), s = sprintf('%d', v);
else, s = sprintf('%.6g', v);
end
end
