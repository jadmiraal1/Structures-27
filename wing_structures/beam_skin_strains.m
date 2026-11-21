function [e_top, e_bot] = beam_skin_strains(P, B, outlines, model, with_drag)
% BEAM_SKIN_STRAINS  Extreme axial skin strains along the beam grid B.s.
%   P        : station properties (station_props / zoned_props)
%   B        : halfspan_beam result for the same P
%   outlines : section loops, one per station of P
%   model    : 'horizontal'  neutral axis taken as horizontal through the
%                            skin centroid: strain = kappa (y - yc) at the
%                            highest and lowest point of the section. Good
%                            when the section's principal axes are within a
%                            few degrees of horizontal.
%              'unsymmetric' neutral axis from unsymmetric bending, the same
%                            assumption the beam stiffness EIxx_eff already
%                            makes (free to deflect sideways). Needed when
%                            the section is inclined or strongly cambered.
%   with_drag: true to add the bending from chordwise loads (B.M_lag).
%              'unsymmetric' only.
% e_top is the strain on the side that is in tension when the wing hogs
% (tip loaded down), e_bot on the other side; + = tension.
hold_i = @(f) interp1(P.s, P.(f), min(max(B.s, P.s(1)), P.s(end)));
if strcmpi(model, 'horizontal')
    yc = hold_i('y_c');
    e_top =  B.kappa .* (hold_i('y_top') - yc);
    e_bot = -B.kappa .* (yc - hold_i('y_bot'));
    return
elseif ~strcmpi(model, 'unsymmetric')
    error('beam_skin_strains:model', 'cfg.strain_model must be ''horizontal'' or ''unsymmetric''.');
end
EIxx = hold_i('EIxx_s') + hold_i('EI_p'); EIxz = hold_i('EIxz_s'); EIzz = hold_i('EIzz_s');
M1 = B.M_hog(:); M2 = zeros(size(M1));
if with_drag, M2 = B.M_lag(:); end
% Each grid point lies between two stations: evaluate the strain field on both
% outlines and interpolate, so the result is smooth along the span.
sg = min(max(B.s(:), P.s(1)), P.s(end));
ns = numel(P.s);
kl = ones(size(sg));
if ns > 1, kl = min(max(floor(interp1(P.s(:), (1:ns)', sg)), 1), ns - 1); end   % station below
ku = min(kl + 1, ns);                                                    % station above
w  = zeros(size(sg));
if ns > 1, w = (sg - P.s(kl)) ./ (P.s(ku) - P.s(kl)); end
e_hi = zeros(size(M1)); e_lo = e_hi;
for k = 1:ns
    i = find(kl == k);
    if ~isempty(i)
        [h, l] = skin_strain(outlines{k}, P.x_c(k), P.y_c(k), EIxx(i), EIxz(i), EIzz(i), M1(i), M2(i));
        e_hi(i) = e_hi(i) + (1 - w(i)) .* h; e_lo(i) = e_lo(i) + (1 - w(i)) .* l;
    end
    i = find(ku == k & ku ~= kl);
    if ~isempty(i)
        [h, l] = skin_strain(outlines{k}, P.x_c(k), P.y_c(k), EIxx(i), EIxz(i), EIzz(i), M1(i), M2(i));
        e_hi(i) = e_hi(i) + w(i) .* h; e_lo(i) = e_lo(i) + w(i) .* l;
    end
end
hog = M1 >= 0;                                        % keep the top/bottom meaning when the moment reverses
e_top = e_hi; e_bot = e_lo;
e_top(~hog) = e_lo(~hog); e_bot(~hog) = e_hi(~hog);
end
