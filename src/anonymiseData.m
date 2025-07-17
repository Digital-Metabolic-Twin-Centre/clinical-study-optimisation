function D_noised = anonymiseData(D, epsilon)
%ANONYMISED Adds Laplace noise to D for differential privacy
%   D         = original m × n count matrix (e.g., IMDs × hospitals)
%   epsilon   = privacy budget (e.g., 1.0 for moderate privacy)
%   D_noised  = anonymised output (rounded, non-negative integers)

    if nargin < 2
        error('Usage: anonymiseD(D, epsilon)');
    end

    sensitivity = 1; % For count queries (max change per individual)
    scale = sensitivity / epsilon;

    % Generate Laplace noise
    noise = laprnd(size(D), 0, scale);

    % Add noise, round to integer, clip negatives
    D_noised = max(round(D + noise), 0);
end

function r = laprnd(sz, mu, b)
%LAPRND Generates Laplace noise (zero-centered, scale b)
%   sz = output size (e.g., size(D)), mu = mean, b = scale

    u = rand(sz) - 0.5;
    r = mu - b * sign(u) .* log(1 - 2 * abs(u));
end
