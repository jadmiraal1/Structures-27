function varargout = wing_config(wing)
% WING_CONFIG  Config of the wing being analysed. Every main script calls
% this instead of a particular <wing>_config, so one set of scripts runs
% any wing.
%
%   wing_config('rear_wing')   make the rear wing the active wing
%   wing_config('front_wing')  back to the front wing (the default)
%   cfg = wing_config()        config of the active wing
%
% The choice is remembered until you change it or restart MATLAB (it
% survives the "clear" at the top of the main scripts). Results files are
% named after the wing (results/<name>_*.mat), so both sets can coexist.
% To add a wing, write <name>_config.m next to front_wing_config.m.
persistent active
if nargin > 0
    if ~exist([wing '_config'], 'file')
        error('wing_config:unknown', 'No %s_config.m found.', wing);
    end
    active = wing;
end
if isempty(active), active = 'front_wing'; end
if nargout > 0 || nargin == 0
    varargout{1} = feval([active '_config']);
else
    fprintf('Active wing: %s\n', active);
end
end
