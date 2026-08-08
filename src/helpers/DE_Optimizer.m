function [xbest, fbest, info] = DE_Optimizer(lb, ub, N_layer, DE, HTLM,...
    Layer_param, initial_param, Target)
% Differential Evolution optimiser (JADE-style current-to-pbest/1 with adaptive
% crossover) that minimises the phase-velocity misfit.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

    % -------------------- Setup & defaults --------------------
    D         = numel(lb);               % decision dimension
    NP        = DE.NP;                 % population size
    CR0       = DE.CR;                 % initial crossover prob (for adaptive mean)
    MaxGen    = DE.Maxitr;             % max generations
    InitPop   = DE.InitPop;            % NP x D initial population
    Verbose   = DE.Verbose;            % print early stop message or not
    TolStall  = DE.TolStall;           % relative stall tolerance
    StallGens = DE.StallGens;          % stop after this many stalled gens
    pbestFrac = DE.pbestFrac;          % fraction of population forming the p-best elite pool

    % -------------------- Initial population ------------------
    pop = max(lb, min(ub, InitPop));     % NP x D

    % -------------------- Initial fitness ---------------------
    fit = zeros(NP,1);                   % fitness for each population member
    parfor i = 1:NP % parfor
        fit(i) = PM_Inv_Objective(pop(i,:), N_layer, HTLM, Layer_param, initial_param, Target);
    end

    % -------------------- Best-so-far bookkeeping --------------
    [fbest, ibest] = min(fit);   % current best value & index
    xbest = pop(ibest,:);        % current best vector

    bestHist = zeros(MaxGen,1);  % best f per generation
    meanHist = bestHist;
    noImprov = 0;                % stall counter
    evals    = NP;               % function eval count (we just evaluated NP)

    % -------------------- Adaptive CR state --------------------
    muCR    = min(1, max(0, CR0));  % running CR mean, always in [0,1]
    sigmaCR = 0.1;                  % sampling std for CRi ~ N(muCR, sigmaCR)
    c_adapt = 0.1;                  % learning rate for updating muCR

    % -------------------- Progress plot init ------------------
    % If a GUI supplied an axes (global WMI_RUN.convAx), draw into it; otherwise use figure 111.
    global WMI_RUN %#ok<GVMIS>
    guiActive = ~isempty(WMI_RUN) && isfield(WMI_RUN,'convAx') && isgraphics(WMI_RUN.convAx);
    if guiActive
        axC = WMI_RUN.convAx; cla(axC,'reset');
    else
        figure(111); clf; axC = gca;
    end
    hLine = plot(axC, NaN, NaN, '.b-', 'LineWidth', 2); hold(axC,'on')
    mLine = plot(axC, NaN, NaN, '.r-', 'LineWidth', 2);
    set(axC, 'YScale', 'log');   % log-scale y-axis
    xlabel(axC,'Generation'); ylabel(axC,'Best misfit'); title(axC,'DE Convergence Progress');
    legend(axC, 'Best', 'Mean');

    % ====================== Main DE loop =======================
    for gen = 1:MaxGen

        % ---- Build the p-best pool once per generation ----
        [~, ord]  = sort(fit, 'ascend');   % sort by fitness (smaller is better)
        pcount    = max(1, round(pbestFrac * NP)); % Defines how many "elite" individuals belong to the p-best pool.

        % ---- Allocate mutant population (one mutant per target) ----
        pop_v = zeros(size(pop));          % will hold v_i for each i

        % ---------------- Mutation: current-to-pbest/1 ----------------
        for i = 1:NP
            % pick an elite individual (p-best) for guidance
            % JADE style: vi​=xi​+F⋅(xpbest​−xi​)+F⋅(xr1​−xr2​); so we need index for xp-best
            pidx = ord( randi(pcount) );             % randomly chooses one solution from the top p-best pool. It’s the “guide” for mutation.
            cands = setdiff(1:NP, [i, pidx]);        % pick two distinct random indices, not i or pidx; 
            % i (the current individual) is excluded because it’s already the base in the formula.
            % pidx (the p-best elite) is excluded because it’s already the guide.
            rr    = cands(randperm(numel(cands), 2));
            r1    = rr(1);  r2 = rr(2);

            Fmax = 0.7; Fmin = 0.3; % Higher F_i: exploration, lower: exploitation
            Fmean = Fmax - (Fmax - Fmin)*(gen/MaxGen);
            F_i = min(Fmax, max(Fmin, Fmean + 0.1*randn));

            % mutant: x_i + F_i*(x_pbest - x_i) + F_i*(x_r1 - x_r2)
            v = pop(i,:) + F_i*(pop(pidx,:) - pop(i,:)) + F_i*(pop(r1,:) - pop(r2,:));
            pop_v(i,:) = v;
        end

        % ------ Crossover + bounds + monotone repair + selection -------
        succCR = nan(NP,1);   % CRi values that succeeded this generation

        parfor i = 1:NP % parfor
            u = pop(i,:);         % trial starts as a copy of parent
            v = pop_v(i,:);       % its mutant
            CRi = muCR + sigmaCR*randn;   % sample from Normal
            if CRi < 0 
                CRi = 0;
            elseif CRi > 1 
                CRi = 1; 
            end
            jrand = randi(D);
            for j = 1:D
                if rand < CRi || j == jrand
                    u(j) = v(j);
                end
            end
            if any(u > ub)
                indx = find(u>ub);
                u(indx) = lb(indx) + rand(size(indx)) .* (ub(indx) - lb(indx));
            end
            if any(u < lb)
                indx = find(u<lb);
                u(indx) = lb(indx) + rand(size(indx)) .* (ub(indx) - lb(indx));
            end

            u = RepairVsProfile(u, lb, ub, Layer_param, initial_param, N_layer); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fu = PM_Inv_Objective(u, N_layer, HTLM, Layer_param, initial_param, Target); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if fu <= fit(i)
                pop(i,:) = u;
                fit(i)   = fu;
                succCR(i) = CRi;   % record CRi that produced an improvement
            end
        end
        evals = evals + NP;  % we just evaluated NP trials

        % ---- Adapt muCR from successful CRi values (if any) ----
        if any(~isnan(succCR))
            muCR = (1 - c_adapt)*muCR + c_adapt * mean(succCR(~isnan(succCR)));
            % keep in [0,1]
            if muCR < 0
                muCR = 0; 
            elseif muCR > 1 
                muCR = 1;
            end
        end

        % ---- Update best-so-far ----
        [fbest, ibest] = min(fit);
        mbest = mean(fit, 'omitnan');
        xbest = pop(ibest,:);

        % ---- Log and plot progress (linear) ----
        bestHist(gen) = fbest;
        meanHist(gen) = mbest;
        set(hLine, 'XData', 1:gen, 'YData', bestHist(1:gen));
        set(mLine, 'XData', 1:gen, 'YData', meanHist(1:gen));
        drawnow limitrate

        % ---- GUI pause / stop control (only when a live GUI axes is present) ----
        if guiActive
            while WMI_RUN.pause && ~WMI_RUN.stop
                drawnow; pause(0.05);
            end
            if WMI_RUN.stop
                if Verbose; fprintf('Stopped by user at gen %d.\n', gen); end
                bestHist = bestHist(1:gen);
                info = struct('history', bestHist, 'gens', gen, 'evals', evals);
                return
            end
        end
        

        % ---- Early stopping on stall (relative improvement test) ----
        fprintf('Gen %4d | best f = %.6g\n', gen, fbest);
        if gen > 1
            % if improvement < TolStall * max(1, |prev best|)
            if abs(bestHist(gen) - bestHist(gen-1)) <= TolStall * max(1, abs(bestHist(gen-1)))
                noImprov = noImprov + 1;
            else
                noImprov = 0;
            end
            if noImprov >= StallGens
                if Verbose
                    fprintf('Early stop at gen %d (stalled %d gens).\n', gen, StallGens);
                end
                bestHist = bestHist(1:gen);
                info = struct('history', bestHist, 'gens', gen, 'evals', evals);
                return
            end
        end

        PM_Inv_PlotProfile(gen, fbest, xbest, N_layer, Layer_param, initial_param, HTLM, Target); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    % -------------------- Outputs --------------------
    info = struct('history', bestHist, 'gens', MaxGen, 'evals', evals);
end
