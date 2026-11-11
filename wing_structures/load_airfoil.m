function af = load_airfoil(source, chord, invert)
% LOAD_AIRFOIL  Airfoil contour scaled to chord, split into surfaces.
%
%   source : path to a Selig-format .dat file (name line, then x y pairs,
%            TE -> upper -> LE -> lower -> TE), OR a 4-char NACA string
%            such as '4412' (placeholder geometry for testing only).
%   chord  : chord length [m]
%   invert : true flips z (wing producing downforce)
%
%   af.xz      : closed contour, N x 2 [m], duplicate TE point removed
%   af.thick_at: function handle, thickness [m] at chordwise x [m]

if numel(source) == 4 && all(isstrprop(source, 'digit'))
    xy = naca4(source, 120);
else
    fid = fopen(source, 'r');
    if fid < 0, error('load_airfoil:file', 'Cannot open %s', source); end
    fgetl(fid);                        % name line
    xy = fscanf(fid, '%f %f', [2 Inf])';
    fclose(fid);
    if any(xy(:) > 1.5)
        error('load_airfoil:format', ['Looks like Lednicer format (point ' ...
              'counts in header). Convert to Selig format.']);
    end
end

if norm(xy(1,:) - xy(end,:)) < 1e-9
    xy(end,:) = [];
end
xy = xy * chord;
if invert, xy(:,2) = -xy(:,2); end

[~, iLE] = min(xy(:,1));
s1 = xy(1:iLE, :);          % TE -> LE (first surface)
s2 = xy(iLE:end, :);        % LE -> TE (second surface)
s1 = flipud(s1);

af.xz = xy;
af.chord = chord;
af.thick_at = @(x) abs(interp1(s1(:,1), s1(:,2), x) - interp1(s2(:,1), s2(:,2), x));
af.surf_at  = @(x) sort([interp1(s1(:,1), s1(:,2), x), interp1(s2(:,1), s2(:,2), x)]);
end

function xy = naca4(code, n)
m = str2double(code(1))/100;
p = str2double(code(2))/10;
tt = str2double(code(3:4))/100;
beta = linspace(0, pi, n)';
x = (1 - cos(beta))/2;
yt = 5*tt*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1036*x.^4);
yc = zeros(size(x)); dyc = zeros(size(x));
if m > 0
    f = x < p;
    yc(f)  = m/p^2*(2*p*x(f) - x(f).^2);
    dyc(f) = 2*m/p^2*(p - x(f));
    yc(~f) = m/(1-p)^2*((1 - 2*p) + 2*p*x(~f) - x(~f).^2);
    dyc(~f)= 2*m/(1-p)^2*(p - x(~f));
end
th = atan(dyc);
xu = x - yt.*sin(th); yu = yc + yt.*cos(th);
xl = x + yt.*sin(th); yl = yc - yt.*cos(th);
xy = [flipud([xu yu]); [xl(2:end) yl(2:end)]];
end
