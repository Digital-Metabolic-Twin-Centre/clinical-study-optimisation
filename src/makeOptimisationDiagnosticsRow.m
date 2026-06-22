function row = makeOptimisationDiagnosticsRow(modelContext, problem, solution, solverName)
% makeOptimisationDiagnosticsRow Build one convergence/feasibility diagnostics row.

    if nargin < 4 || isempty(solverName)
        solverName = "not returned";
    end

    x = getSolutionVector(solution);

    solverStatus = getFirstReturnedString(solution, ...
        ["status", "stat", "origStat", "solverStatus", "exitflag"]);
    finalObjectiveValue = getFirstReturnedNumber(solution, ...
        ["obj", "objective", "objectiveValue", "f", "full"]);
    maximumFeasibilityResidual = computeMaximumFeasibilityResidual(problem, x);

    row = table( ...
        string(modelContext), ...
        string(solverName), ...
        solverStatus, ...
        finalObjectiveValue, ...
        maximumFeasibilityResidual, ...
        "N/A", ...
        'VariableNames', { ...
            'model_context', ...
            'solver_name', ...
            'solver_status', ...
            'final_objective_value', ...
            'maximum_feasibility_residual', ...
            'mip_gap'});
end

function x = getSolutionVector(solution)
    if isstruct(solution) && isfield(solution, 'xyz') && ~isempty(solution.xyz)
        x = full(solution.xyz(:));
    else
        x = NaN;
    end
end

function value = getFirstReturnedString(solution, fieldNames)
    value = "not returned";
    if ~isstruct(solution)
        return
    end
    for i = 1:numel(fieldNames)
        fieldName = char(fieldNames(i));
        if isfield(solution, fieldName) && ~isempty(solution.(fieldName))
            value = string(solution.(fieldName));
            return
        end
    end
end

function value = getFirstReturnedNumber(solution, fieldNames)
    value = NaN;
    if ~isstruct(solution)
        return
    end
    for i = 1:numel(fieldNames)
        fieldName = char(fieldNames(i));
        if isfield(solution, fieldName) && ~isempty(solution.(fieldName))
            candidate = solution.(fieldName);
            if isnumeric(candidate) && isscalar(candidate)
                value = double(candidate);
                return
            end
        end
    end
end

function residual = computeMaximumFeasibilityResidual(problem, x)
    if numel(x) ~= numel(problem.c)
        residual = NaN;
        return
    end

    residuals = 0;
    AxMinusB = full(problem.A * x - problem.b(:));
    csense = upper(cellstr(problem.csense(:)));

    eqMask = strcmp(csense, 'E');
    if any(eqMask)
        residuals = max(residuals, max(abs(AxMinusB(eqMask))));
    end

    leMask = strcmp(csense, 'L');
    if any(leMask)
        residuals = max(residuals, max(max(AxMinusB(leMask), 0)));
    end

    geMask = strcmp(csense, 'G');
    if any(geMask)
        residuals = max(residuals, max(max(-AxMinusB(geMask), 0)));
    end

    lowerViolation = full(problem.lb(:)) - x;
    upperViolation = x - full(problem.ub(:));
    residuals = max(residuals, max([lowerViolation; upperViolation; 0]));

    residual = residuals;
end
