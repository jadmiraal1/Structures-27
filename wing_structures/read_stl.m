function T = read_stl(fname, scale)
% READ_STL  Read a binary or ASCII STL into an nT x 9 triangle array.
%   T(i,:) = [x1 y1 z1 x2 y2 z2 x3 y3 z3], multiplied by SCALE
%   (e.g. 1e-3 for an STL exported in millimetres).
if nargin < 2, scale = 1; end

fid = fopen(fname, 'r');
if fid < 0, error('read_stl:file', 'Cannot open %s', fname); end
fseek(fid, 0, 'eof'); nbytes = ftell(fid); fseek(fid, 0, 'bof');
hdr = fread(fid, 80, 'uint8=>char')';
ntri = fread(fid, 1, 'uint32');

if ~isempty(ntri) && nbytes == 84 + 50*ntri
    % Binary: per triangle 12 float32 (normal + 3 vertices) + uint16
    raw = fread(fid, [25, ntri], 'uint16=>uint16');   % 50 bytes = 25 uint16
    fclose(fid);
    raw = raw(1:24, :);                                % drop attribute word
    f = typecast(raw(:), 'single');
    f = reshape(f, 12, ntri)';
    T = double(f(:, 4:12));
else
    fclose(fid);
    txt = fileread(fname);
    tok = regexp(txt, 'vertex\s+(\S+)\s+(\S+)\s+(\S+)', 'tokens');
    if isempty(tok)
        error('read_stl:format', 'Not a recognised binary or ASCII STL: %s (header "%s")', fname, strtrim(hdr));
    end
    v = str2double(vertcat(tok{:}));
    T = reshape(v', 9, [])';
end
T = T * scale;
end
