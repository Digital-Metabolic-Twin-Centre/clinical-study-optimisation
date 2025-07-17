function [U_heur, H_heur, costT_heur, costG_heur, totalCost_heur] = ...
         HeuristicBaselineAllocation(y, T, G, r, a, B, capT, capG)
% HeuristicBaselineAllocation with optional sample caps per platform.
% Ensures T*U ≥ 2S and G*H ≥ 2S, and optionally sum(U,2) ≤ capT and sum(H,2) ≤ capG.
%
% Inputs:
%   y     (m×1) – patients per disease
%   T     (p×q) – targeted platform biomarker coverage (0/1)
%   G     (p×w) – global   platform biomarker coverage (0/1)
%   r     (q×1) – targeted cost per sample
%   a     (w×1) – global   cost per sample
%   B     (p×m) – biomarker-disease incidence (0/1)
%   capT  (q×1) – (optional) max samples per targeted platform
%   capG  (w×1) – (optional) max samples per global platform
%
% Outputs:
%   U_heur, H_heur – sample allocation matrices
%   costT_heur, costG_heur, totalCost_heur – cost breakdown

    y = y(:);  % ensure column
    [p, m] = size(B);
    [~, q] = size(T);
    [~, w] = size(G);
    S2 = 2 * (B * diag(y));  % double coverage requirement

    % flip
    r = r';
    a = a';

    if nargin < 7 || isempty(capT), capT = inf(q,1); end
    if nargin < 8 || isempty(capG), capG = inf(w,1); end

    %% Allocate targeted samples
    U_heur = zeros(q, m);
    usedT  = zeros(q,1);  % track total usage
    for d = 1:m
        unmet = S2(:,d) - T * U_heur(:,d);
        while any(unmet > 0)
            coverage = sum(T .* (unmet > 0), 1)';  % 1×q
            feasible = (coverage > 0) & (usedT < capT);  % respect cap
            if ~any(feasible)
                error('No feasible targeted platform left for disease %d', d);
            end
            cost_eff = inf(1, q);
            cost_eff(feasible) = r(feasible) ./ coverage(feasible);
            [~, jbest] = min(cost_eff);
            inds = find(T(:, jbest) & (unmet > 0));
            k = ceil(max(unmet(inds)));
            max_allow = capT(jbest) - usedT(jbest);
            k = min(k, max_allow);
            if k <= 0
                feasible(jbest) = false;
                continue;
            end
            U_heur(jbest, d) = U_heur(jbest, d) + k;
            usedT(jbest)     = usedT(jbest) + k;
            unmet = S2(:,d) - T * U_heur(:,d);
        end
    end

    %% Allocate global samples
    H_heur = zeros(w, m);
    usedG  = zeros(w,1);
    for d = 1:m
        unmet = S2(:,d) - G * H_heur(:,d);
        while any(unmet > 0)
            coverage = sum(G .* (unmet > 0), 1)';  % 1×w
            feasible = (coverage > 0) & (usedG < capG);
            if ~any(feasible)
                error('No feasible global platform left for disease %d', d);
            end
            cost_eff = inf(1, w);
            cost_eff(feasible) = a(feasible) ./ coverage(feasible);
            [~, kbest] = min(cost_eff);
            inds = find(G(:, kbest) & (unmet > 0));
            k = ceil(max(unmet(inds)));
            max_allow = capG(kbest) - usedG(kbest);
            k = min(k, max_allow);
            if k <= 0
                feasible(kbest) = false;
                continue;
            end
            H_heur(kbest, d) = H_heur(kbest, d) + k;
            usedG(kbest)     = usedG(kbest) + k;
            unmet = S2(:,d) - G * H_heur(:,d);
        end
    end

    %% Compute costs
    costT_heur     = sum((sum(U_heur, 2)) .* r);
    costG_heur     = sum((sum(H_heur, 2)) .* a);
    totalCost_heur = costT_heur + costG_heur;

    %% Coverage checks
    assert(all(all(T * U_heur >= S2)), 'Targeted coverage assertion failed');
    assert(all(all(G * H_heur >= S2)), 'Global coverage assertion failed');
    assert(all(sum(U_heur,2) <= capT + 1e-6), 'Targeted platform cap violated');
    assert(all(sum(H_heur,2) <= capG + 1e-6), 'Global platform cap violated');
end
