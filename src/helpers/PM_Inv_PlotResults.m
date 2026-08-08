function PM_Inv_PlotResults(files, path, opts, axDisp, axProf)
% Plots saved PM_Inv inversion results (dispersion curves + Vs-depth profiles)
% from the selected .mat files into the given dispersion / profile axes.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026
%
% files : cell array of .mat filenames (or a single char)   path : their folder
% axDisp/axProf : target axes for dispersion / profile (a 1x2 subplot figure is
%                 created if they are omitted).
%
% opts fields (all optional):
%   Dmax   - maximum plotting depth (m)                         [30]
%   Nbest  - models drawn per trial when Mode = 'All'           [1]
%   Mode   - selection mode: 'All' | 'Best' | 'MeanBest'        ['All']
%   nModes - number of theoretical modes to overlay (1..5)      [3]
%
% Mode meanings:
%   'All'      : the Nbest lowest-misfit models from EACH selected trial file.
%   'Best'     : the single lowest-misfit model per h_count (across its trials).
%   'MeanBest' : the mean of the best-per-trial models for each h_count.
%
% Layout-robust: model rows are read as h=3, Vs=4, modes=end-1, predominant=end,
% so it works for both the viscoelastic (9-row param) and elastic (8-row) results.

    if ischar(files); files = {files}; end
    if nargin < 3 || isempty(opts); opts = struct(); end
    if nargin < 5 || ~isgraphics(axDisp) || ~isgraphics(axProf)
        figure('Color','w','Name','PM_Inv Post-Processing');
        axDisp = subplot(1,2,1); axProf = subplot(1,2,2);
    end
    Dmax   = getdef(opts,'Dmax',30);
    Nbest  = getdef(opts,'Nbest',1);
    Mode   = getdef(opts,'Mode','All');
    nModes = getdef(opts,'nModes',3);

    modeColors = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; ...
                  0.494 0.184 0.556; 0.301 0.745 0.933];

    % ---------------- load + index every file ----------------
    R = struct('hc',{},'tr',{},'param',{},'Target',{},'mis',{},'mR',{},'dR',{},'hr',{});
    for k = 1:numel(files)
        try
            d = load(fullfile(path, files{k}));
        catch ME
            warning('Could not load %s (%s)', files{k}, ME.message); continue;
        end
        if ~isfield(d,'param') || ~isfield(d,'Target')
            warning('Skipping %s: missing param/Target.', files{k}); continue;
        end
        tok = regexp(files{k}, 'h_count_([-0-9.eE]+)_trial_([-0-9.eE]+)', 'tokens', 'once');
        if isempty(tok); hc = k; tr = 1; else; hc = str2double(tok{1}); tr = str2double(tok{2}); end
        nR   = size(d.param,1);
        good = ~cellfun(@isempty, d.param(2,:));       % keep populated generations
        hr   = '';
        if isfield(d,'Layer_param') && isfield(d.Layer_param,'h_range'); hr = d.Layer_param.h_range; end
        r = struct('hc',hc, 'tr',tr, 'param',{d.param(:,good)}, 'Target',d.Target, ...
                   'mis',cell2mat(d.param(2,good)), 'mR',nR-1, 'dR',nR, 'hr',hr);
        R(end+1) = r; %#ok<AGROW>
    end
    if isempty(R); warning('No valid PM_Inv result files selected.'); return; end

    Tg   = R(1).Target.fv_file;
    freq = Tg(:,1);
    hcs  = [R.hc];
    uH   = unique(hcs, 'stable');
    grpColors = lines(max(numel(uH),1));

    % ---------------- assemble the list of models to draw ----------------
    % Each P(i): zStep,vsStep (profile), modes (Nf x m), domi (Nf x 1), hc
    P = struct('zStep',{},'vsStep',{},'modes',{},'domi',{},'hc',{},'misfits',{});
    switch lower(Mode)
        case 'best'
            for u = 1:numel(uH)
                ri = find(hcs==uH(u));  bestVal = inf; bR = ri(1); bG = 1;
                for r = ri
                    [mv, gi] = min(R(r).mis);
                    if mv < bestVal; bestVal = mv; bR = r; bG = gi; end
                end
                P(end+1) = makeEntry(R(bR), bG, Dmax); %#ok<AGROW>
            end

        case 'meanbest'
            zq = linspace(0, Dmax, 400);
            for u = 1:numel(uH)
                ri = find(hcs==uH(u));
                VS = []; DOMI = []; MOD = []; MISF = [];
                for r = ri
                    [~, gi] = min(R(r).mis);
                    m = makeEntry(R(r), gi, Dmax);
                    VS(end+1,:)    = stairsVs(m.zStep, m.vsStep, zq);      %#ok<AGROW>
                    DOMI(end+1,:)  = padvec(m.domi(:).', numel(freq));     %#ok<AGROW>
                    MOD(:,:,end+1) = padmodes(m.modes, nModes, numel(freq)); %#ok<AGROW>
                    MISF(end+1)    = m.misfits;                            %#ok<AGROW>
                end
                e = struct('zStep',zq, 'vsStep',mean(VS,1,'omitnan'), ...
                           'modes',mean(MOD,3,'omitnan'), 'domi',mean(DOMI,1,'omitnan').', ...
                           'hc',uH(u), 'misfits',MISF);
                P(end+1) = e; %#ok<AGROW>
            end

        otherwise   % 'All'
            for r = 1:numel(R)
                for g = topGens(R(r).mis, Nbest)
                    P(end+1) = makeEntry(R(r), g, Dmax); %#ok<AGROW>
                end
            end
    end

    % ---------------- draw ----------------
    cla(axDisp,'reset'); cla(axProf,'reset');
    hold(axDisp,'on'); hold(axProf,'on');

    profLegH = gobjects(1,numel(uH));   % one profile handle per h_count (legend)
    for i = 1:numel(P)
        gi  = find(uH == P(i).hc, 1);
        col = grpColors(gi,:);

        % --- profile (Vs vs depth) ---
        hp = plot(axProf, P(i).vsStep, P(i).zStep, '-', 'Color', col, 'LineWidth', 1.3);
        if ~isgraphics(profLegH(gi)); profLegH(gi) = hp; end

        % --- dispersion: theoretical modes (thin) then predominant (thick) ---
        km = min(nModes, size(P(i).modes,2));
        for mm = 1:km
            plot(axDisp, freq, P(i).modes(:,mm), '-', 'Color', modeColors(mm,:), 'LineWidth', 0.5);
        end
        plot(axDisp, freq, P(i).domi, '*', 'Color', col, 'LineWidth', 1);
    end

    % target dispersion (with std error bars if present)
    if size(Tg,2) >= 3
        hT = errorbar(axDisp, Tg(:,1), Tg(:,2), Tg(:,3), '.k', 'MarkerSize', 8, 'CapSize', 3);
    else
        hT = plot(axDisp, Tg(:,1), Tg(:,2), '.k', 'MarkerSize', 8);
    end

    % ---------------- cosmetics ----------------
    axis(axProf,'ij'); box(axProf,'on'); ylim(axProf,[0 Dmax]);
    xlabel(axProf,'Shear-wave velocity (m/s)'); ylabel(axProf,'Depth (m)');
    set(axProf,'FontName','times','FontSize',11,'TickDir','out','LineWidth',1);
    title(axProf,'Vs profiles');
    valid = isgraphics(profLegH);
    if any(valid)
        pref = 'h';
        if strcmpi(R(1).hr,'LR_based'); pref = 'LR'; elseif strcmpi(R(1).hr,'LN_based'); pref = 'LN'; end
        vi = find(valid); labs = cell(1, numel(vi));
        for q = 1:numel(vi)
            gi = vi(q);
            ms = [P([P.hc] == uH(gi)).misfits];     % misfits of the displayed models for this h_count
            if isempty(ms)
                labs{q} = sprintf('%s = %g', pref, uH(gi));
            elseif max(ms) - min(ms) < 1e-9
                labs{q} = sprintf('%s = %g   misfit %.3g', pref, uH(gi), min(ms));
            else
                labs{q} = sprintf('%s = %g   misfit [%.3g, %.3g]', pref, uH(gi), min(ms), max(ms));
            end
        end
        legend(axProf, profLegH(valid), labs, 'Location','best');
    end

    box(axDisp,'on'); set(axDisp,'XScale','log','FontName','times','FontSize',11,'TickDir','out','LineWidth',1);
    xlabel(axDisp,'Frequency (Hz)'); ylabel(axDisp,'Phase velocity (m/s)');
    xlim(axDisp,[min(freq) max(freq)]); title(axDisp,'Dispersion');
    % legend proxies: modes + predominant + target
    proxy = gobjects(0); lab = {};
    for mm = 1:min(nModes,size(modeColors,1))
        proxy(end+1) = plot(axDisp, NaN, NaN, '-', 'Color', modeColors(mm,:), 'LineWidth', 1.5); %#ok<AGROW>
        lab{end+1}   = sprintf('Mode %d', mm); %#ok<AGROW>
    end
    proxy(end+1) = plot(axDisp, NaN, NaN, '*', 'Color', [.4 .4 .4], 'LineWidth', 1.5); lab{end+1} = 'Predominant';
    proxy(end+1) = hT; lab{end+1} = 'Target';
    legend(axDisp, proxy, lab, 'Location','best');
end

% ===================== local helpers =====================
function v = getdef(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function gens = topGens(mis, N)
% generations of the N smallest DISTINCT misfits
    [mu, ia] = unique(mis(:).');           % sorted ascending, first index of each
    N = min(N, numel(mu));
    gens = ia(1:N).';
end

function P = makeEntry(R, g, Dmax)
    h  = R.param{3, g};   vs = R.param{4, g};
    [zStep, vsStep] = makeStep(h, vs, Dmax);
    P = struct('zStep',zStep, 'vsStep',vsStep, ...
               'modes',R.param{R.mR, g}, 'domi',R.param{R.dR, g}, 'hc',R.hc, 'misfits',R.mis(g));
end

function [zStep, vsStep] = makeStep(h, vs, Dmax)
    h = h(:).'; vs = vs(:).';
    h(isnan(h)) = [];
    dpth   = cumsum(h);
    vsStep = repelem(vs, 2);
    zStep  = [0 repelem(dpth,2) Dmax];
    n = min(numel(vsStep), numel(zStep));   % guard against length mismatch
    zStep = zStep(1:n); vsStep = vsStep(1:n);
end

function vq = stairsVs(zStep, vsStep, zq)
% resample a step profile (zStep,vsStep) onto query depths zq (piecewise-constant).
% zStep contains duplicate depths (a staircase), so interp1 can't be used; take
% the Vs of the last step whose depth is <= the query depth.
    vq = zeros(size(zq));
    for j = 1:numel(zq)
        idx = find(zStep <= zq(j), 1, 'last');
        if isempty(idx); idx = 1; end
        vq(j) = vsStep(idx);
    end
end

function v = padvec(v, n)
    v = v(:).';
    if numel(v) < n; v(end+1:n) = NaN; else; v = v(1:n); end
end

function M = padmodes(modes, k, nf)
% return nf x k matrix (first k modes; pad with NaN if fewer / shorter)
    M = nan(nf, k);
    if isempty(modes); return; end
    kk = min(k, size(modes,2));
    rr = min(nf, size(modes,1));
    M(1:rr, 1:kk) = modes(1:rr, 1:kk);
end
