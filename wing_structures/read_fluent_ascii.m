function F = read_fluent_ascii(fname)
% READ_FLUENT_ASCII  Read a Fluent "Export > Solution Data > ASCII" file of
% wall-face data (cell-center location, comma delimited).
%
% Required columns (Fluent names): x/y/z-coordinate, pressure,
% x/y/z-wall-shear, x/y/z-face-area.
%
% Returns struct with fields (one row per face, SI units as exported):
%   xyz [N x 3]  face centroid
%   p   [N x 1]  static (gauge) pressure
%   tau [N x 3]  wall shear stress vector
%   A   [N x 3]  face area vector
%   Amag[N x 1]  face area

fid = fopen(fname, 'r');
if fid < 0, error('read_fluent_ascii:file', 'Cannot open %s', fname); end
hdr = fgetl(fid);
fclose(fid);
names = strtrim(strsplit(hdr, ','));
names = lower(names);

need = {'x-coordinate','y-coordinate','z-coordinate','pressure', ...
        'x-wall-shear','y-wall-shear','z-wall-shear', ...
        'x-face-area','y-face-area','z-face-area'};
col = zeros(1, numel(need));
for k = 1:numel(need)
    j = find(strcmp(names, need{k}), 1);
    if isempty(j)
        error('read_fluent_ascii:column', ['Column "%s" not found. Header was:\n%s\n' ...
              'Export Static Pressure (gauge), X/Y/Z Wall Shear Stress and X/Y/Z Face Area.'], need{k}, hdr);
    end
    col(k) = j;
end

D = dlmread(fname, ',', 1, 0);
F.xyz  = D(:, col(1:3));
F.p    = D(:, col(4));
F.tau  = D(:, col(5:7));
F.A    = D(:, col(8:10));
F.Amag = sqrt(sum(F.A.^2, 2));
F.file = fname;
end
