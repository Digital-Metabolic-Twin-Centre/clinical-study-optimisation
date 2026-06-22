function row = makeOptimisationDiagnosticsRow(modelContext, problem, solution, solvingTimeSeconds, solverName)
% makeOptimisationDiagnosticsRow Build one convergence/feasibility diagnostics row.

    if nargin < 5 || isempty(solverName)
        solverName = "not returned";
    end

    x = getSolutionVector(solution);

    numberDecisionVariables = numel(problem.c);
    numberEqualityConstraints = countConstraints(problem, 'E');
    numberInequalityConstraints = countConstraints(problem, {'G', 'L'});
    numberBoundConstraints = countFiniteBounds(problem);
    numberObjectiveVariables = nnz(logical(problem.p)) + nnz(logical(problem.q));
    degreesOfFreedom = numberDecisionVariables - numberEqualityConstraints;

    solverStatus = getFirstReturnedString(solution, ...
        ["status", "stat", "origStat", "solverStatus", "exitflag"]);
    stoppingCriterion = getFirstReturnedString(solution, ...
        ["stoppingCriterion", "stopCriterion", "termination", "message", "msg"]);
    numberCardinalityIterations = getFirstReturnedNumber(solution, ...
        ["number_cardinality_iterations", "cardinalityIterations", "iterations", "iter", "itn", "nIter"]);
    finalObjectiveValue = getFirstReturnedNumber(solution, ...
        ["obj", "objective", "objectiveValue", "f", "full"]);
    finalRelativeObjectiveChange = getFinalRelativeObjectiveChange(solution);
    maximumFeasibilityResidual = computeMaximumFeasibilityResidual(problem, x);

    row = table( ...
        string(modelContext), ...
        numberDecisionVariables, ...
        numberEqualityConstraints, ...
        numberInequalityConstraints, ...
        numberBoundConstraints, ...
        numberObjectiveVariables, ...
        degreesOfFreedom, ...
        solvingTimeSeconds, ...
        string(solverName), ...
        solverStatus, ...
        stoppingCriterion, ...
        numberCardinalityIterations, ...
        finalObjectiveValue, ...
        finalRelativeObjectiveChange, ...
        maximumFeasibilityResidual, ...
        "N/A", ...
        'VariableNames', { ...
            'model_context', ...
            'number_decision_variables', ...
            'number_equality_constraints', ...
            'number_inequality_constraints', ...
            'number_bound_constraints', ...
            'number_objective_variables', ...
            'degrees_of_freedom', ...
            'solving_time_seconds', ...
            'solver_name', ...
            'solver_status', ...
            'stopping_criterion', ...
            'number_cardinality_iterations', ...
            'final_objective_value', ...
            'final_relative_objective_change', ...
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

function n = countConstraints(problem, senses)
    if ~isfield(problem, 'csense') || isempty(problem.csense)
        n = NaN;
        return
    end
    csense = upper(cellstr(problem.csense(:)));
    if ischar(senses)
        senses = {senses};
    end
    senses = upper(string(senses));
    n = sum(ismember(string(csense), senses));
end

function n = countFiniteBounds(problem)
    n = 0;
    if isfield(problem, 'lb') && ~isempty(problem.lb)
        n = n + nnz(isfinite(full(problem.lb(:))));
    end
    if isfield(problem, 'ub') && ~isempty(problem.ub)
        n = n + nnz(isfinite(full(problem.ub(:))));
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

function value = getFinalRelativeObjectiveChange(solution)
    value = NaN;
    history = getObjectiveHistory(solution);
    if numel(history) < 2
        return
    end
    previousObjective = history(end - 1);
    finalObjective = history(end);
    value = abs(finalObjective - previousObjective) ./ max(1, abs(previousObjective));
end

function history = getObjectiveHistory(solution)
    history = [];
    if ~isstruct(solution)
        return
    end

    fieldNames = ["objectiveHistory", "objHistory", "objective_history", "history"];
    for i = 1:numel(fieldNames)
        fieldName = char(fieldNames(i));
        if ~isfield(solution, fieldName) || isempty(solution.(fieldName))
            continue
        end

        candidate = solution.(fieldName);
        if isnumeric(candidate) && isvector(candidate)
            history = double(candidate(:));
            return
        end

        if isstruct(candidate)
            nested = getFirstNumericVectorField(candidate, ["obj", "objective", "objectiveValue"]);
            if ~isempty(nested)
                history = nested;
                return
            end
        end
    end
end

function values = getFirstNumericVectorField(s, fieldNames)
    values = [];
    for i = 1:numel(fieldNames)
        fieldName = char(fieldNames(i));
        if isfield(s, fieldName) && isnumeric(s.(fieldName)) && isvector(s.(fieldName))
            values = double(s.(fieldName)(:));
            return
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
