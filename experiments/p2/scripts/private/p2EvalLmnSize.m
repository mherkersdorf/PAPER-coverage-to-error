function M = p2EvalLmnSize(N)
%P2EVALLMNSIZE  The number of local models at dataset size N: M = max(1, round(min(2.25 sqrt(N), N/10))).
%
%   Fixed per size, identical for every
%   family. 201 at N = 8000; the cap N/10 binds below N of about 530.
%
%   Inputs:
%       N  (1 x n) double  dataset sizes
%   Outputs:
%       M  (1 x n) double  local models
M = max(1, round(min(2.25 * sqrt(N), N / 10)));
end
