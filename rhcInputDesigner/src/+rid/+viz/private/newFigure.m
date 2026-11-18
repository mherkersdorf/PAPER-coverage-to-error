function fig = newFigure(name, visible)
%NEWFIGURE  One figure, named after the run, visible only when asked.
%
%   Every +viz entry point creates its own figure rather than drawing into
%   `gcf`: a plotting call that silently overwrites whatever the user had open
%   is not reproducible, and a test cannot assert against a figure it does not
%   own.
%
%   Inputs:
%       name     (1x1) string                        the figure's Name
%       visible  (1x1) matlab.lang.OnOffSwitchState
%   Outputs:
%       fig  (1x1) matlab.ui.Figure
fig = figure("Name", name, "NumberTitle", "off", "Visible", visible, ...
             "Color", "w");

end
