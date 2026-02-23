function [U_greedy, H_greedy, costT, costG, totalCost, stats] = ...
    GlobalGreedyMarginalSampleAllocation(y, T, G, r, a, B, capT, capG, opts)
% GlobalGreedyMarginalSampleAllocation (ENHANCED + RANDOM RESTARTS)
% Capacity-aware global greedy benchmark for sample allocation.
%
% This version adds RANDOM RESTARTS to make batchMode=false robust:
% - If batchMode=false (1-sample steps), greedy can get stuck even when a feasible
%   solution exists. Random restarts (random tie-breaking among maxima) reduce
%   this failure mode while keeping the *same model/algorithm family*.
%
% Inputs: (same as before) plus optional in opts:
%   opts.nRestarts   (default 1)  : number of random restarts (only used when batchMode=false)
%   opts.rngSeed     (default 1)  : base RNG seed for reproducibility across restarts
%   opts.pickBestBy  (default 'totalCost') : 'totalCost' or 'targetedThenGlobal'
%
% Outputs:
%   U_greedy (q×m), H_greedy (w×m)
%   costT, costG, totalCost
%   stats struct with diagnostics

    if nargin < 9 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'batchMode'),          opts.batchMode = true; end
    if ~isfield(opts, 'maxIter'),            opts.maxIter = 2e7; end
    if ~isfield(opts, 'verbose'),            opts.verbose = false; end
    if ~isfield(opts, 'activationPenaltyT'), opts.activationPenaltyT = 0; end
    if ~isfield(opts, 'activationPenaltyG'), opts.activationPenaltyG = 0; end
    if ~isfield(opts, 'preferFewPlatforms'), opts.preferFewPlatforms = true; end
    if ~isfield(opts, 'tieBreak'),           opts.tieBreak = 'minIndex'; end
    if ~isfield(opts, 'useTwoXCoverage'),    opts.useTwoXCoverage = true; end

    % --- Random restart options ---
    if ~isfield(opts, 'nRestarts'),          opts.nRestarts = 1; end
    if ~isfield(opts, 'rngSeed'),            opts.rngSeed = 1; end
    if ~isfield(opts, 'pickBestBy'),         opts.pickBestBy = 'totalCost'; end

    y = y(:);
    [p, m] = size(B);
    [pT, q] = size(T);
    [pG, w] = size(G);
    assert(pT == p && pG == p, 'T and G must have same biomarker dimension as B.');

    if nargin < 7 || isempty(capT), capT = inf(q,1); end
    if nargin < 8 || isempty(capG), capG = inf(w,1); end
    capT = capT(:); capG = capG(:);

    r = r(:);
    a = a(:);

    % Coverage requirement matrix (p×m)
    if opts.useTwoXCoverage
        Sreq = 2 * (B * diag(y));
    else
        Sreq = (B * diag(y));
    end

    % Decide how many restarts to run:
    % - If batchMode=true, the greedy is already less path-dependent; restarts are usually unnecessary.
    % - If batchMode=false, restarts help avoid "stuck" trajectories.
    if opts.batchMode
        nRuns = 1;
    else
        nRuns = max(1, round(opts.nRestarts));
    end

    best = struct();
    best.found = false;
    best.totalCost = inf;
    best.costT = inf;
    best.costG = inf;
    best.U = [];
    best.H = [];
    best.statsRun = struct();

    % ===========================
    %   Random-restart wrapper
    % ===========================
    for run = 1:nRuns
        if ~opts.batchMode
            rng(opts.rngSeed + run - 1, 'twister'); % reproducible restarts
        end

        % Each run may use random tie-breaking among maxima when batchMode=false.
        % We keep opts.tieBreak as-is if batchMode=true.
        optsRun = opts;
        if ~opts.batchMode
            optsRun.tieBreak = 'random';
        end

        % Run targeted + global greedies; if either fails, this run fails.
        try
            [U_run, statsT_run] = solveOneSideGreedy( ...
                'targeted', Sreq, T, r, capT, optsRun);

            [H_run, statsG_run] = solveOneSideGreedy( ...
                'global',   Sreq, G, a, capG, optsRun);

        catch ME
            if opts.verbose
                fprintf('[Restart %d/%d] failed: %s\n', run, nRuns, ME.message);
            end
            continue;
        end

        % Costs
        costT_run = sum(sum(U_run, 2) .* r);
        costG_run = sum(sum(H_run, 2) .* a);
        totalCost_run = costT_run + costG_run;

        % Keep best
        take = false;
        if ~best.found
            take = true;
        else
            switch lower(opts.pickBestBy)
                case 'totalcost'
                    take = totalCost_run < best.totalCost;
                case 'targetedthenglobal'
                    if costT_run < best.costT
                        take = true;
                    elseif costT_run == best.costT
                        take = costG_run < best.costG;
                    end
                otherwise
                    take = totalCost_run < best.totalCost;
            end
        end

        if take
            best.found = true;
            best.totalCost = totalCost_run;
            best.costT = costT_run;
            best.costG = costG_run;
            best.U = U_run;
            best.H = H_run;

            best.statsRun = struct();
            best.statsRun.restartUsed = run;
            best.statsRun.targeted = statsT_run;
            best.statsRun.global   = statsG_run;
        end
    end

    if ~best.found
        error('GlobalGreedyMarginalAllocation: no feasible solution found after %d run(s).', nRuns);
    end

    % Outputs
    U_greedy = best.U;
    H_greedy = best.H;
    costT = best.costT;
    costG = best.costG;
    totalCost = best.totalCost;

    % Assertions (feasibility)
    assert(all(all(T * U_greedy >= Sreq)), 'Targeted coverage assertion failed');
    assert(all(all(G * H_greedy >= Sreq)), 'Global coverage assertion failed');
    assert(all(sum(U_greedy,2) <= capT + 1e-6), 'Targeted platform cap violated');
    assert(all(sum(H_greedy,2) <= capG + 1e-6), 'Global platform cap violated');

    % Stats
    stats = best.statsRun;
    stats.nRunsTried = nRuns;
    stats.bestTotalCost = best.totalCost;
    stats.bestCostT = best.costT;
    stats.bestCostG = best.costG;
    stats.requirementMatrixUsed = Sreq;

end

% ========================================================================
% Local helper: solve one side (Targeted or Global) with same greedy logic
% ========================================================================
function [Alloc, statsSide] = solveOneSideGreedy(sideName, Sreq, M, c, cap, opts)
% sideName: 'targeted' or 'global' (for messages only)
% M: coverage matrix (p×k), k = q or w
% c: costs (k×1)
% cap: caps (k×1)

    [p, m] = size(Sreq);
    [pM, k] = size(M);
    assert(pM == p, '%s: coverage matrix biomarker dimension mismatch.', sideName);

    Alloc = zeros(k, m);
    remCap = cap(:);
    usedPlatform = false(k,1);

    R = max(Sreq - M * Alloc, 0); % unmet (p×m)

    iter = 0;
    while any(R(:) > 0)
        iter = iter + 1;
        if iter > opts.maxIter
            error('%s greedy exceeded maxIter (%g).', sideName, opts.maxIter);
        end

        feasPlat = (remCap > 0);
        if ~any(feasPlat)
            error('No feasible %s platform left (all caps exhausted).', sideName);
        end

        gainAmount = (M' * R);          % k×m
        gainAmount(~feasPlat, :) = 0;

        if all(gainAmount(:) <= 0)
            error('%s greedy: no remaining feasible gain (coverage gap or caps too tight).', ...
                capitalize(sideName));
        end

        score = gainAmount ./ c(:);     % k×m

        % Activation penalty (optional)
        if strcmpi(sideName, 'targeted') && opts.activationPenaltyT > 0
            if opts.preferFewPlatforms
                penaltyVec = opts.activationPenaltyT * double(~usedPlatform);
            else
                penaltyVec = opts.activationPenaltyT * ones(k,1);
            end
            score = score - penaltyVec * ones(1,m);

        elseif strcmpi(sideName, 'global') && opts.activationPenaltyG > 0
            if opts.preferFewPlatforms
                penaltyVec = opts.activationPenaltyG * double(~usedPlatform);
            else
                penaltyVec = opts.activationPenaltyG * ones(k,1);
            end
            score = score - penaltyVec * ones(1,m);
        end

        score(~feasPlat, :) = -inf;

        % Pick best (platform, disease)
        [bestVal, linIdx] = max(score(:));
        if ~isfinite(bestVal)
            error('%s greedy: no remaining feasible move (all scores -Inf).', capitalize(sideName));
        end

        if strcmpi(opts.tieBreak, 'random')
            maxMask = (score(:) == bestVal);
            idxs = find(maxMask);
            linIdx = idxs(randi(numel(idxs)));
        end

        [pBest, dBest] = ind2sub([k, m], linIdx);

        % Step size
        if opts.batchMode
            coveredIdx = (M(:, pBest) > 0) & (R(:, dBest) > 0);
            if ~any(coveredIdx)
                step = 1;
            else
                step = ceil(max(R(coveredIdx, dBest)));
            end
        else
            step = 1;
        end

        step = min(step, remCap(pBest));
        if step <= 0
            remCap(pBest) = 0;
            continue;
        end

        Alloc(pBest, dBest) = Alloc(pBest, dBest) + step;
        remCap(pBest)       = remCap(pBest) - step;
        usedPlatform(pBest) = true;

        % Update unmet only for this disease
        R(:, dBest) = max(Sreq(:, dBest) - M * Alloc(:, dBest), 0);

        if opts.verbose && mod(iter, 5000) == 0
            fprintf('[%s] iter=%d, remaining-unmet=%g, activePlatforms=%d\n', ...
                capitalize(sideName), iter, nnz(R > 0), nnz(usedPlatform));
        end
    end

    % return stats
    statsSide = struct();
    statsSide.iter = iter;
    statsSide.activePlatforms = find(usedPlatform);
    statsSide.nActivePlatforms = nnz(usedPlatform);
    statsSide.remainingCap = remCap;
end

function s = capitalize(s)
    if isempty(s), return; end
    s(1) = upper(s(1));
end