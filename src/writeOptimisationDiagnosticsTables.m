function writeOptimisationDiagnosticsTables(diagnosticsTable, resultFolder)
% writeOptimisationDiagnosticsTables Write convergence diagnostics as CSV and LaTeX.

    if ~exist(resultFolder, 'dir')
        mkdir(resultFolder);
    end

    csvPath = fullfile(resultFolder, 'optimisation_convergence_diagnostics.csv');
    texPath = fullfile(resultFolder, 'optimisation_convergence_diagnostics.tex');

    writetable(diagnosticsTable, csvPath);
    writeDiagnosticsLatexTable(diagnosticsTable, texPath);
end

function writeDiagnosticsLatexTable(T, texPath)
    fid = fopen(texPath, 'w');
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '\\begin{tabular}{%s}\n', repmat('l', 1, width(T)));
    fprintf(fid, '\\hline\n');
    fprintf(fid, '%s \\\\\n', strjoin(escapeLatex(T.Properties.VariableNames), ' & '));
    fprintf(fid, '\\hline\n');

    for i = 1:height(T)
        values = strings(1, width(T));
        for j = 1:width(T)
            values(j) = formatTableValue(T{i, j});
        end
        fprintf(fid, '%s \\\\\n', strjoin(escapeLatex(values), ' & '));
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
end

function value = formatTableValue(value)
    if iscell(value)
        value = value{1};
    end

    if isstring(value) || ischar(value)
        value = string(value);
        return
    end

    if isnumeric(value)
        if isnan(value)
            value = "NaN";
        elseif isscalar(value) && abs(value - round(value)) < eps(max(1, abs(value)))
            value = string(sprintf('%.0f', value));
        else
            value = string(sprintf('%.6g', value));
        end
        return
    end

    value = string(value);
end

function escaped = escapeLatex(values)
    escaped = string(values);
    replacements = [
        "\", "\textbackslash{}";
        "_", "\_";
        "%", "\%";
        "&", "\&";
        "#", "\#";
        "{", "\{";
        "}", "\}";
        "$", "\$"
    ];

    for i = 1:size(replacements, 1)
        escaped = replace(escaped, replacements(i, 1), replacements(i, 2));
    end
end
