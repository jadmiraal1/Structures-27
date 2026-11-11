% MAIN_STAGE1  Section properties: skin shell vs spar options.
%
% Question this stage answers: how stiff is the skin shell compared with
% each spar candidate, and which spar is lightest for a given stiffness?
%
% Because the skin is bonded to the ribs only (not along the spar), the ribs
% enforce matching vertical deflection and twist at rib stations but do not
% transfer spanwise shear. So the skin and spar each bend about their OWN
% neutral axis, and the combined stiffness is approximately
%     EI_total ~ EI_skin + EI_spar      (not a parallel-axis composite)
% Stage 2 (coupled two-beam FE) checks this and gives the rib transfer loads.
%
% Units: SI. All inputs marked PLACEHOLDER need your values.
clear; close all;

%% ---------------- INPUTS ----------------
M = materials_input();

geo.airfoil  = '4412';     % PLACEHOLDER: path to Selig .dat, e.g. 'mainplane.dat'
geo.chord    = 0.30;       % PLACEHOLDER [m]
geo.invert   = true;       % downforce
geo.spar_xc  = 0.30;       % PLACEHOLDER spar centre location, x/c
geo.gap      = 0.5e-3;     % bond-line / fit clearance each side [m]

skin_layup   = [0 45 0];   % PLACEHOLDER twill skin layup (deg from spanwise)

% Rectangular tube stock: [width_chordwise  depth  wall] in inches. PLACEHOLDER list.
tube_in = [0.50 0.75 0.049;
           0.50 0.75 0.065;
           0.50 1.00 0.049;
           0.50 1.00 0.065;
           0.75 1.00 0.065;
           0.50 1.25 0.065;
           1.00 1.00 0.065];

% Sandwich web options
core_t      = [6.35e-3 9.525e-3];               % 1/4", 3/8" core
face_layups = {[0], [45], [45 45], [0 45 0]};   % per face
cap_plies   = 0:4;                               % 0-deg twill plies per cap

%% ---------------- SETUP ----------------
af = load_airfoil(geo.airfoil, geo.chord, geo.invert);
if numel(geo.airfoil) == 4
    fprintf('*** PLACEHOLDER GEOMETRY: NACA %s ***\n', geo.airfoil);
end
x_spar = geo.spar_xc * geo.chord;

bound_name = {'LOW stiffness bound', 'HIGH stiffness bound'};
results = cell(1, 2);

for k = 1:2
    tw = bound_props(M.twill, k);
    al = bound_props(M.alum,  k);

    %% Skin shell
    % Airfoil outer mould line used as the skin midline. Error is O(t/thickness),
    % small for a ~1 mm skin in a ~35 mm deep section.
    skin = laminate_clt(tw, skin_layup);
    Ssk  = thinwall_closed_section(af.xz, skin.h, skin.Ex, skin.Gxy, skin.rho);

    % Usable depth inside the skin, minimum over a chordwise width w
    depth_avail = @(w) min(af.thick_at(linspace(x_spar - w/2, x_spar + w/2, 25))) ...
                       - 2*skin.h - 2*geo.gap;

    R = struct('name', {}, 'family', {}, 'mass', {}, 'EI', {}, 'GJ', {}, 'fits', {});

    %% Aluminum tubes
    for i = 1:size(tube_in, 1)
        b = tube_in(i,1)*0.0254; h = tube_in(i,2)*0.0254; t = tube_in(i,3)*0.0254;
        S = thinwall_closed_section(rect_tube_midline(b, h, t), t, al.E, al.G, al.rho);
        R(end+1) = struct('name', sprintf('Al %.2fx%.2fx%.3f in', tube_in(i,:)), ...
            'family', 1, 'mass', S.mass_per_len, 'EI', S.EIxx, 'GJ', S.GJ, ...
            'fits', h <= depth_avail(b)); %#ok<SAGROW>
    end

    %% Sandwich webs (full available depth)
    for ic = 1:numel(core_t)
        for jf = 1:numel(face_layups)
            face = laminate_clt(tw, face_layups{jf});
            tpanel = core_t(ic) + 2*face.h;
            htot = depth_avail(tpanel);
            for nc = cap_plies
                if nc == 0
                    cap = []; fam = 2;
                else
                    cap = struct('lam', laminate_clt(tw, zeros(1, nc)), 'w', tpanel);
                    fam = 3;
                end
                S = sandwich_web_section(htot, face, core_t(ic), M.nomex.rho, cap, M.bond_areal);
                R(end+1) = struct('name', sprintf('SW core %.1fmm face [%s] caps %d', ...
                    core_t(ic)*1e3, num2str(face_layups{jf}), nc), 'family', fam, ...
                    'mass', S.mass_per_len, 'EI', S.EIxx, 'GJ', S.GJ, 'fits', true); %#ok<SAGROW>
            end
        end
    end

    %% Report
    fprintf('\n================ %s ================\n', bound_name{k});
    fprintf('Skin [%s]: t = %.2f mm, Ex = %.1f GPa, Gxy = %.1f GPa\n', ...
        num2str(skin_layup), skin.h*1e3, skin.Ex/1e9, skin.Gxy/1e9);
    fprintf('  EIxx = %.0f N m^2   EIxx_eff (camber-coupled) = %.0f N m^2\n', Ssk.EIxx, Ssk.EIxx_eff);
    fprintf('  GJ   = %.0f N m^2   mass = %.3f kg/m\n', Ssk.GJ, Ssk.mass_per_len);
    fprintf('  Usable spar depth at x/c = %.2f: %.1f mm\n', geo.spar_xc, depth_avail(0)*1e3);
    fprintf('\n%-38s %8s %10s %9s %9s %6s\n', 'Spar option', 'kg/m', 'EI [Nm2]', 'EI/EIsk', 'GJ/GJsk', 'fits');
    [~, order] = sort([R.mass]); yn = {'NO', 'yes'};
    for i = order
        fprintf('%-38s %8.3f %10.0f %9.3f %9.4f %6s\n', R(i).name, R(i).mass, R(i).EI, ...
            R(i).EI/Ssk.EIxx_eff, R(i).GJ/Ssk.GJ, yn{R(i).fits + 1});
    end
    results{k} = struct('skin', Ssk, 'R', R);
end

%% Plot: mass per length vs bending stiffness (LOW bound)
R = results{1}.R; Ssk = results{1}.skin;
figure; hold on; grid on;
mk = {'s', 'o', '^'}; lbl = {'Al rect tube', 'Sandwich web', 'Sandwich web + twill caps'};
for fam = 1:3
    sel = [R.family] == fam & [R.fits];
    if any(sel), plot([R(sel).EI], [R(sel).mass], mk{fam}, 'DisplayName', lbl{fam}); end
end
yl = ylim;
plot(Ssk.EIxx_eff*[1 1], yl, 'k--', 'DisplayName', 'Skin shell EI');
xlabel('Spar EI_{xx} [N m^2]'); ylabel('Spar mass [kg/m]');
title('Spar options, low stiffness bound (spars that fit only)');
legend('Location', 'northwest');
print('-dpng', 'stage1_spar_tradeoff.png');
