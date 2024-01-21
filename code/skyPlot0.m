function hpol = skyPlot0(varargin)
%% initial 
red_my =[0.8500 0.3250 0.0980];
bule_my = [0 0.4470 0.7410];
yellow_my = [0.9290 0.6940 0.1250];
purple_my =[0.4940 0.1840 0.5560];
el_cutoff = 15;

%% Check arguments and sort them ==========================================
[hAxis, args, nargs] = axescheck(varargin{:});

if nargs < 3 || nargs > 5
    error('Requires 3 or 5 data arguments.')
elseif nargs == 3
    [az, el, prn]   = deal(args{1:3});
    cn0 = [];
    flag = true(length(prn), 1);
    line_style      = 'auto';    
elseif nargs == 4
    [az, el, prn, cn0] = deal(args{1:4});
    flag = true(length(prn), 1);
    line_style      = 'auto';
elseif nargs == 5
    [az, el, prn, cn0, flag] = deal(args{1:4});
    line_style      = 'auto';    
end

if ischar(az) || ischar(el) || ischar(prn)
    error('AZ and EL must be numeric.');
end

if ~isequal(size(az), size(el))
    error('AZ and EL must be same size.');
end

%% Prepare axis ===========================================================
hAxis = newplot(hAxis);

%--- Get x-axis text color so grid is in same color -----------------------
tc = get(hAxis, 'xcolor');

hold(hAxis, 'on');

%--- Plot white background ------------------------------------------------
rectangle('position', [-90, -90, 180, 180], ...
          'Curvature', [1 1], ...
          'facecolor', 'white', ...
          'edgecolor', tc);
%% Plot spokes ============================================================

%--- Find spoke angles ----------------------------------------------------
% Only 6 lines are needed to divide circle into 12 parts
th = (1:6) * 2*pi / 12;

%--- Convert spoke end point coordinate to Cartesian system ---------------
cst = cos(th); snt = sin(th);
cs = [cst; -cst];
sn = [snt; -snt];

%--- Plot the spoke lines -------------------------------------------------
line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ...
    'handlevisibility', 'off');

%% Annotate spokes in degrees =============================================
rt = 1.1 * 90;

for i = 1:max(size(th))

    %--- Write text in the first half of the plot -------------------------
    text(rt*snt(i), rt*cst(i), int2str(i*30), ...
        'horizontalalignment', 'center', 'handlevisibility', 'off');

    if i == max(size(th))
        loc = int2str(0);
    else
        loc = int2str(180 + i*30);
    end

    %--- Write text in the opposite half of the plot ----------------------
    text(-rt*snt(i), -rt*cst(i), loc, ...
        'handlevisibility', 'off', 'horizontalalignment', 'center');
end

%% Plot elevation grid ====================================================

%--- Define a "unit" radius circle ----------------------------------------
th = 0 : pi/50 : 2*pi;
xunit = cos(th);
yunit = sin(th);

%--- Plot elevation grid lines and tick text ------------------------------
for elevation = 0 : 15 : 90
    elevationSpherical = 90*cos((pi/180) * elevation);

    if elevation == el_cutoff
        line(yunit * elevationSpherical, xunit * elevationSpherical, ...
            'lineStyle', ':', 'color', red_my, 'linewidth', 1.5, ...
            'handlevisibility', 'off');
    else
        line(yunit * elevationSpherical, xunit * elevationSpherical, ...
        'lineStyle', ':', 'color', tc, 'linewidth', 1, ...
        'handlevisibility', 'off');
    end

    text(0, elevationSpherical, num2str(elevation), ...
        'BackgroundColor', 'none', 'horizontalalignment','center', ...
        'handlevisibility', 'off');
end

%--- Set view to 2-D ------------------------------------------------------
view(0, 90);

%--- Set axis limits ------------------------------------------------------
%save some space for the title
axis([-95 95 -90 101]);

%% Transform elevation angle to a distance to the center of the plot ------
elSpherical = 90*cos(el * pi/180);

%--- Transform data to Cartesian coordinates ------------------------------
yy = elSpherical .* cos(az * pi/180);
xx = elSpherical .* sin(az * pi/180);

%% Plot data on top of the grid ===========================================

if strcmp(line_style, 'auto')
    %--- Plot with "default" line style -----------------------------------
%      hpol = plot(hAxis, xx', yy', '.-');
    %hpol = plot(hAxis, xx', yy', line_style);
else
    %--- Plot with user specified line style ------------------------------
    % The same line style and color will be used for all satellites
    hpol = plot(hAxis, xx', yy', line_style);
end


%--- Mark the position of the satellite ------------------------------
plot(hAxis, xx(flag), yy(flag), 'p', 'color', red_my, 'MarkerSize', 10 , 'Markerfacecolor' ,red_my);
plot(hAxis, xx(~flag), yy(~flag), 'p', 'color', red_my, 'MarkerSize', 10 , 'Markerfacecolor' ,yellow_my);
if nargs == 5
    legend('use in PVT', 'only tracked');
    legend('boxoff');
    legend('Location','northeastoutside');
end

%--- Place satellite PRN numbers and CN0 -------------------
for i = 1:length(prn)
    text(xx(i), yy(i), num2str(prn(i), 'C%02d'), 'color', bule_my);
    if (~isempty(cn0))
        text(xx(i), yy(i)-4, int2str(cn0(i)), 'color', purple_my);
    end
end

%--- Make sure both axis have the same data aspect ratio ------------------
axis(hAxis, 'equal');

%--- Switch off the standard Cartesian axis -------------------------------
axis(hAxis, 'off');
end