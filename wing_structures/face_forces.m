function [Ff, Fp, Fv] = face_forces(F)
% FACE_FORCES  Force on each CFD wall face [N], in the CFD frame.
%   pressure part  Fp = p * A      (A = exported face area vector)
%   viscous part   Fv = tau * |A|
% Sign convention verified against Fluent Report > Forces: summing these
% over the front-wing zone reproduces Fluent's pressure and viscous
% totals (see main_loads check).
Fp = F.p .* F.A;
Fv = F.tau .* F.Amag;
Ff = Fp + Fv;
end
