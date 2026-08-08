function figPos = GetFigurePosition(widthScale, heightScale, leftDiv, bottomDiv)
% Computes a screen-relative [left bottom width height] position vector so figures
% are sized and placed consistently across different screen resolutions.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026


    if nargin < 4
        error('Usage: GetFigurePosition(widthScale, heightScale, leftDiv, bottomDiv)');
    end

    screen_size = get(0, 'ScreenSize');
    fig_width   = screen_size(3) * widthScale;
    fig_height  = screen_size(4) * heightScale;
    fig_left    = (screen_size(3) - fig_width) / leftDiv;
    fig_bottom  = (screen_size(4) - fig_height) / bottomDiv;

    figPos = [fig_left, fig_bottom, fig_width, fig_height];
end