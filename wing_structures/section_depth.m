function S = section_depth(xy, nq)
% SECTION_DEPTH  Chord and thickness distribution of a closed airfoil loop.
%   xy : M x 2 loop [streamwise, vertical]. Streamwise axis points FORWARD,
%        so the leading edge is the max-x point and the trailing edge the
%        min-x point.
%   nq : number of chordwise query points (default 101)
%
% Returns chord, LE/TE points, chord-line angle, and thickness t(xc) at
% xc = 0..1 measured along the streamwise projection from the LE.
% Thickness is the vertical distance between the two surfaces at a given
% streamwise station (the depth available to a vertical spar web).
if nargin < 2, nq = 101; end
[~, iLE] = max(xy(:,1));
[~, iTE] = min(xy(:,1));
M = size(xy, 1);
% two chains from LE to TE, going each way round the loop
idx1 = mod((iLE-1):(iLE-1+mod(iTE-iLE, M)), M) + 1;
idx2 = mod((iLE-1):-1:(iLE-1-mod(iLE-iTE, M)), M) + 1;
c1 = xy(idx1, :); c2 = xy(idx2, :);

xLE = xy(iLE,1); xTE = xy(iTE,1);
xq = linspace(xLE, xTE, nq)';
y1 = interp_chain(c1, xq);
y2 = interp_chain(c2, xq);

S.chord_x = xLE - xTE;                      % streamwise projected chord
S.LE = xy(iLE,:); S.TE = xy(iTE,:);
S.chord = norm(S.LE - S.TE);
S.angle_deg = atan2d(S.LE(2) - S.TE(2), S.LE(1) - S.TE(1));  % + = nose up
S.xc = (xLE - xq) / (xLE - xTE);
S.x  = xq;
S.y_hi = max(y1, y2);
S.y_lo = min(y1, y2);
S.t  = S.y_hi - S.y_lo;
[S.t_max, im] = max(S.t);
S.xc_tmax = S.xc(im);
end

function yq = interp_chain(c, xq)
% chain runs LE -> TE; keep the first occurrence of each x so a chain that
% doubles back near the LE stays single valued
[xs, ia] = unique(c(:,1), 'stable');
ys = c(ia, 2);
[xs, o] = sort(xs); ys = ys(o);
yq = interp1(xs, ys, xq, 'linear', 'extrap');
end
