function [Y, y_byDisease, totalPatients, stats] = GlobalGreedyMarginalPatientRecruitment(D_k, patientNumberUB, HCPcap_k, opts)
% GlobalGreedyMarginalPatientRecruitment
% One-by-one global greedy recruitment across diseases and hospitals.
% At each step choose (disease i, hospital j) that maximizes:
%   score(i,j) = scarcity(i) * remainingAvailability(i,j)
% subject to:
%   - Y(i,j) < D_k(i,j)
%   - hospital remaining capacity > 0
%   - disease remaining demand > 0
%
% Inputs:
%   D_k            (m×n_k) availability matrix for selected hospitals
%   patientNumberUB (m×1) max patients per disease
%   HCPcap_k       (n_k×1) max patients per hospital
%   opts (struct, optional)
%       .useScarcity  true (default)
%       .maxSteps     default sum(patientNumberUB)
%       .tieBreak     'minIndex' (default) or 'random'
%       .rngSeed      default 1 (used if tieBreak='random')
%       .verbose      false (default)
%
% Outputs:
%   Y             (m×n_k) allocation
%   y_byDisease   (m×1) recruited patients per disease
%   totalPatients scalar total recruited
%   stats struct with diagnostics

    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'useScarcity'), opts.useScarcity = true; end
    if ~isfield(opts,'maxSteps'),    opts.maxSteps    = []; end
    if ~isfield(opts,'tieBreak'),    opts.tieBreak    = 'minIndex'; end
    if ~isfield(opts,'rngSeed'),     opts.rngSeed     = 1; end
    if ~isfield(opts,'verbose'),     opts.verbose     = false; end

    patientNumberUB = patientNumberUB(:);
    HCPcap_k = HCPcap_k(:);

    [m, n_k] = size(D_k);
    assert(numel(patientNumberUB) == m, 'patientNumberUB length must match rows of D_k.');
    assert(numel(HCPcap_k) == n_k, 'HCPcap_k length must match cols of D_k.');

    if isempty(opts.maxSteps)
        opts.maxSteps = sum(patientNumberUB);
    end

    if strcmpi(opts.tieBreak, 'random')
        rng(opts.rngSeed, 'twister');
    end

    Y = zeros(m, n_k);
    hospitalRem = HCPcap_k(:)';          % 1×n_k
    diseaseRem  = patientNumberUB(:);    % m×1

    totalAvailPerDisease = sum(D_k, 2);  % m×1
    if opts.useScarcity
        scarcity = 1 ./ max(1, totalAvailPerDisease);
    else
        scarcity = ones(m,1);
    end

    step = 0;
    while step < opts.maxSteps
        step = step + 1;

        needDisease = (diseaseRem > 0);  % m×1

        feas = (D_k > Y) & (hospitalRem > 0) & needDisease; % m×n_k
        if ~any(feas(:)), break; end

        remAvail = (D_k - Y); % m×n_k

        score = (scarcity * ones(1, n_k)) .* remAvail;
        score(~feas) = -inf;

        bestVal = max(score(:));
        if ~isfinite(bestVal), break; end

        linIdx = find(score(:) == bestVal);

        if numel(linIdx) > 1 && strcmpi(opts.tieBreak, 'random')
            linIdx = linIdx(randi(numel(linIdx)));
        else
            linIdx = linIdx(1); % deterministic: first max in linear order
        end

        [iBest, jBest] = ind2sub([m, n_k], linIdx);

        Y(iBest, jBest) = Y(iBest, jBest) + 1;
        hospitalRem(jBest) = hospitalRem(jBest) - 1;
        diseaseRem(iBest)  = diseaseRem(iBest)  - 1;

        if opts.verbose && mod(step, 5000) == 0
            fprintf('[GreedyRecruit] step=%d, total=%d\n', step, sum(Y(:)));
        end
    end

    y_byDisease   = sum(Y, 2);
    totalPatients = sum(y_byDisease);

    stats = struct();
    stats.steps = step;
    stats.remainingHospitalCapacity = hospitalRem(:);
    stats.remainingDiseaseDemand    = diseaseRem(:);
end