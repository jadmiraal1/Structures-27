function K = kink_corners(s, y_top, y_bot, N_top, N_bot, curv_min)
% KINK_CORNERS  Find where the skin turns in front view (loft corners) and
% the sideways load that the turn puts on the skin.
%   s            : spanwise positions [m], uniform spacing
%   y_top, y_bot : height of the top and bottom skin lines in front view [m]
%   N_top, N_bot : spanwise skin force per unit chord width [N/m], + tension
%   curv_min     : |d(theta)/ds| above which the skin counts as turning [rad/m]
%
% A skin carrying force N that turns through an angle dtheta pushes
% sideways with N*dtheta per unit chord width (like a belt on a pulley).
% Each contiguous turning region is one corner. Returns, per corner:
%   s       : position of the sharpest turn [m]
%   s_range : [start end] of the turning region [m]
%   dtheta  : angle change of the top line [rad]
%   F_top, F_bot : |N| * |dtheta| at the corner [N/m of chord]
s = s(:); ds = s(2) - s(1);
sm = @(v) conv(v(:), ones(5,1)/5, 'same');
yt = sm(y_top); yb = sm(y_bot);
yt(1:2) = y_top(1:2); yt(end-1:end) = y_top(end-1:end);
yb(1:2) = y_bot(1:2); yb(end-1:end) = y_bot(end-1:end);
tht = atan(gradient(yt, ds)); thb = atan(gradient(yb, ds));
kt = gradient(tht, ds); kb = gradient(thb, ds);
turn = max(abs(kt), abs(kb)) > curv_min;

d = diff([0; turn; 0]);
i0 = find(d == 1); i1 = find(d == -1) - 1;
K = struct('s', {}, 's_range', {}, 'dtheta', {}, 'F_top', {}, 'F_bot', {});
for k = 1:numel(i0)
    a = max(i0(k) - 1, 1); b = min(i1(k) + 1, numel(s));
    [~, j] = max(max(abs(kt(a:b)), abs(kb(a:b)))); j = a + j - 1;
    dtt = tht(b) - tht(a); dtb = thb(b) - thb(a);
    K(end+1) = struct('s', s(j), 's_range', [s(a) s(b)], 'dtheta', dtt, ...
        'F_top', abs(N_top(j)) * abs(dtt), 'F_bot', abs(N_bot(j)) * abs(dtb)); %#ok<AGROW>
end
end
