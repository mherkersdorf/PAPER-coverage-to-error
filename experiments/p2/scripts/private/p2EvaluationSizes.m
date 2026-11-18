function sizes = p2EvaluationSizes()
%P2EVALUATIONSIZES  The dataset sizes of A and B: 80, evenly spaced on a log axis from 50 to 8000.
%
%   factor ~1.066 between neighbours, rounded to
%   integers, distinct. The arms and the APRBS are evaluated on prefixes of their
%   one 8000-sample run; the multisine and the chirp get one signal per size.
%   Used by the generation and the evaluation, so both use
%   ONE list.
%
%   Outputs:
%       sizes  (1 x 80) double  increasing integers, sizes(1) = 50, sizes(end) = 8000
sizes = round(logspace(log10(50), log10(8000), 80));
if numel(unique(sizes)) ~= 80 || sizes(1) ~= 50 || sizes(end) ~= 8000
    error("experiments:p2EvaluationSizes:notDistinct", ...
        "The rounded sizes are not 80 distinct integers from 50 to 8000.");
end
end
