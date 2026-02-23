function [Y, y_byDisease, totalPatients, stats] = HeuristicPatientRecruitment(D_k, patientNumberUB, HCPcap_k, opts)
% HeuristicPatientRecruitment
% Per-disease proportional allocation + greedy fill, respecting:
%   - disease upper bounds patientNumberUB(i)
%   - hospital capacity caps HCPcap_k(j)
%   - hospital availability D_k(i,j)
%
% Inputs:
%   D_k            (m×n_k) availability matrix for selected hospitals
%   patientNumberUB (m×1) max patients to recruit per disease
%   HCPcap_k       (n_k×1) max patients per hospital
%   opts (struct, optional)
%       .rounding  'round' (default) or 'floor' or 'ceil'
%       .verbose   false (default)
%
% Outputs:
%   Y             (m×n_k) allocation
%   y_byDisease   (m×1) recruited patients per disease
%   totalPatients scalar total recruited
%   stats struct with diagnostics

    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'rounding'), opts.rounding = 'round'; end
    if ~isfield(opts,'verbose'),  opts.verbose  = false; end

    patientNumberUB = patientNumberUB(:);
    HCPcap_k = HCPcap_k(:);

    [m, n_k] = size(D_k);
    assert(numel(patientNumberUB) == m, 'patientNumberUB length must match rows of D_k.');
    assert(numel(HCPcap_k) == n_k, 'HCPcap_k length must match cols of D_k.');

    Y = zeros(m, n_k);
    hospitalRem = HCPcap_k(:)'; % 1×n_k remaining capacity

    for i = 1:m
        diseasePreval = D_k(i,:);
        totalAvail = sum(diseasePreval);
        if totalAvail <= 0 || patientNumberUB(i) <= 0
            continue;
        end

        wts = diseasePreval ./ max(eps, totalAvail);
        raw = wts * patientNumberUB(i);

        switch lower(opts.rounding)
            case 'round'
                rawAlloc = round(raw);
            case 'floor'
                rawAlloc = floor(raw);
            case 'ceil'
                rawAlloc = ceil(raw);
            otherwise
                error('Unknown rounding option: %s', opts.rounding);
        end

        feasible = min(rawAlloc, diseasePreval);
        feasible = min(feasible, hospitalRem);

        assigned  = sum(feasible);
        remaining = patientNumberUB(i) - assigned;

        % Greedy fill: add one-by-one where still possible
        while remaining > 0
            candidates = (diseasePreval > feasible) & (hospitalRem > feasible);
            if ~any(candidates), break; end

            canImprove = (diseasePreval - feasible) .* candidates;
            [bestVal, bestH] = max(canImprove);
            if bestVal <= 0, break; end

            feasible(bestH)   = feasible(bestH) + 1;
            hospitalRem(bestH)= hospitalRem(bestH) - 1;
            remaining         = remaining - 1;
        end

        Y(i,:) = feasible;
        hospitalRem = hospitalRem - feasible;

        if opts.verbose
            fprintf('[Heuristic] disease %d recruited %d/%d\n', i, sum(feasible), patientNumberUB(i));
        end
    end

    y_byDisease   = sum(Y, 2);
    totalPatients = sum(y_byDisease);

    stats = struct();
    stats.remainingHospitalCapacity = hospitalRem(:);
    stats.usedHospitalCapacity      = HCPcap_k - stats.remainingHospitalCapacity;
end