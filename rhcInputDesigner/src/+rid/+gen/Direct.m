classdef Direct < rid.gen.Generator
%DIRECT  Every future input sample is a decision variable.
%
%   The default parametrization and the reference against which the deferred
%   reduced ones (B-spline, step) are judged. There is no candidate SET: the
%   decision space is `nHorizon * nInputs` dimensional, so this class hands the
%   optimizer a DESCRIPTION of the feasible set and, later, a handful of starting
%   points -- never a list.
%
%   THE DECISION VECTOR IS CHANNEL-MAJOR:
%
%       theta = uSeq(:)          uSeq is (nHorizon x nInputs)
%       theta((c-1)*nHorizon + j) is channel c at horizon step j
%
%   so each channel owns a CONTIGUOUS block of `nHorizon` entries. This ordering is
%   not free: the rate constraint is `kron(eye(nInputs), D)` with `D` the
%   first-difference matrix, and a
%   block-diagonal `kron` acts on contiguous per-channel blocks. A sample-major
%   ordering would need `kron(D, eye(nInputs))` and a different `b`, so the two
%   choices are not interchangeable and the next reader should not have to re-derive
%   which one section 3 assumes.
%
%   `nHorizon` IS HELD HERE AND NOWHERE ELSE. It is a required constructor argument
%   with no default anywhere, on the precedent: the method description makes it user-defined and fixed for the run; the time constant is
%   guidance and never a formula, so no defensible default exists and inventing one
%   would hide the decision. The generator is the SINGLE HOLDER -- the loop
%   reads `gen.nHorizon` rather than keeping its own copy, and the surrogate's
%   simulation length, the representation's future coordinates and the penalty's
%   `gamma` all read this one value.
%
%   The channel count is likewise not stored: it is read from the constraints, which
%   read it from the representation.
%
%   Construction:
%       gen = rid.gen.Direct(cons, 20)
%
%   Inputs:
%       constraints  (1x1) rid.gen.Constraints  the declared admissible region
%       nHorizon     (1x1) double               prediction horizon, samples
%
%   PHYSICAL UNITS throughout: +gen works in physical units and normalization
%   happens at the representation boundary.
%
%   THE CONTRACT IT IMPLEMENTS IS `rid.gen.Generator`, which is where the
%   interface every parametrization satisfies is written down. This class is the
%   default one and, until reduced parametrizations exist, the only one.
%
%   THE GRADIENT COMPOSITION HELPERS ARE NOT HERE. `rid.gen.horizonGradient` and
%   `rid.gen.gradientToTheta` were statics of this class and are now package
%   functions: selecting horizon rows and converting a physical gradient
%   through a diagonal factor is generic composition maths and belongs to no
%   parametrization. `horizonWeights` and `targetSequence` stay, because `gamma`
%   and the target ramp ARE this parametrization's.

    properties (SetAccess = immutable)
        % BOTH are declared without validators because both implement abstract
        % properties of rid.gen.Generator, which MATLAB does not permit
        % validators on -- the same constraint rid.rep.Narx carries for nDim.
        % Shapes:
        %
        %   constraints  (1x1) rid.gen.Constraints  the hard box, PHYSICAL
        %   nHorizon     (1x1) double               prediction horizon in samples
        %
        % `constraints` joined the contract at the constraints contract; `nHorizon` has its SINGLE
        % HOLDER rule from the nHorizon ownership rule. The constructor's `arguments` block types both
        % on the way in, which is where the validation went.
        constraints
        nHorizon
    end

    methods

        function obj = Direct(constraints, nHorizon)
            arguments
                constraints (1,1) rid.gen.Constraints
                nHorizon    (1,1) double {mustBeInteger, mustBePositive}
            end
            % nHorizon = 1 is legal and is the MYOPIC BASELINE: greedy one-step
            % design, with no interior rate constraint at all (the band matrix is
            % empty) and only the rate limit against u(k-1) acting. Keeping it a
            % legal configuration rather than an error makes "how much does the
            % horizon buy?" an experiment rather than a code change.
            obj.constraints = constraints;
            obj.nHorizon    = nHorizon;
        end

        function p = nInputs(obj)
            %NINPUTS  Input channels, read from the constraints.
            %
            %   A method and not a property: `p` has one source, the representation,
            %   and a stored copy here could disagree with it.
            %
            %   Outputs:
            %       p  (1x1) double
            p = numel(obj.constraints.channelName);
        end

        function n = nParams(obj)
            %NPARAMS  Length of the decision vector.
            %
            %   Takes no argument: `nHorizon` lives on this object, and a parameter
            %   here would be a second source of one number that could disagree with
            %   the first.
            %
            %   Outputs:
            %       n  (1x1) double  nHorizon * nInputs
            n = obj.nHorizon * obj.nInputs();
        end

        function uSeq = toTrajectory(obj, theta)
            %TOTRAJECTORY  The decision vector as a physical input trajectory.
            %
            %   Inputs:
            %       theta  (nParams x 1) double  channel-major, see the class header
            %   Outputs:
            %       uSeq   (nHorizon x nInputs) double  physical units
            arguments
                obj
                theta (:,1) double {mustBeReal}
            end
            if numel(theta) ~= obj.nParams()
                error("rid:gen:invalidTheta", ...
                    "theta has %d entries; this parametrization has %d " + ...
                    "(nHorizon = %d, nInputs = %d).", numel(theta), ...
                    obj.nParams(), obj.nHorizon, obj.nInputs());
            end
            uSeq = reshape(theta, obj.nHorizon, obj.nInputs());
        end

        function uApp = applied(obj, uSeq)
            %APPLIED  The unit committed this iteration: ONE SAMPLE.
            %
            %   The applied unit is defined by the parametrization, not by the loop.
            %   For a spline parametrization it would be the first segment; here it
            %   is the first sample, and +rhc commits exactly what this returns.
            %
            %   Inputs:
            %       uSeq  (nHorizon x nInputs) double
            %   Outputs:
            %       uApp  (1 x nInputs) double
            arguments
                obj
                uSeq (:,:) double {mustBeReal}
            end
            obj.assertTrajectorySize(uSeq);
            uApp = uSeq(1, :);
        end

        function theta = shift(obj, thetaPrev)
            %SHIFT  The warm start: shift by one, repeat the last sample.
            %
            %   FEASIBLE BY CONSTRUCTION, given that the applied input equals the
            %   previous solution's first sample. Every interior rate constraint
            %   carries over unchanged, the new first variable's constraint against
            %   u(k-1) is the previous solution's own first interior constraint, and
            %   repeating the last sample adds a zero increment -- which the rate box
            %   contains, because rid.gen.Constraints rejects one that does not. No
            %   feasibility repair is ever needed, and none is done here.
            %
            %   A SOFT output constraint can still be violated by the warm start
            %   after a surrogate update, because the trajectory moves although
            %   theta did not. Box and rate constraints are on u alone and are immune.
            %
            %   Inputs:
            %       thetaPrev  (nParams x 1) double
            %   Outputs:
            %       theta      (nParams x 1) double
            arguments
                obj
                thetaPrev (:,1) double {mustBeReal}
            end
            uPrev = obj.toTrajectory(thetaPrev);
            uNext = [uPrev(2:end, :); uPrev(end, :)];
            theta = uNext(:);
        end

        function [lb, ub, A, b] = feasibleSet(obj, hist)
            %FEASIBLESET  The hard constraints, exactly.
            %
            %   STATE-DEPENDENT through u(k-1): the first variable's box is the
            %   amplitude box intersected with `u(k-1) + rate`, because
            %   `|u_hat(k) - u(k-1)| <= duMax` is a BOUND, not a row of A.
            %
            %   The interior rate constraints are the first-difference matrix,
            %   blocked per channel:
            %
            %       D = spdiags([-1 +1], 0:1, nHorizon-1, nHorizon)
            %       A = [ kron(I, D) ; -kron(I, D) ],  b = [ duMax ; -duMin ]
            %
            %   Two blocks and not one because the rate box need not be symmetric:
            %   an actuator that opens faster than it closes is the normal case.
            %
            %   Inputs:
            %       hist  (1x1) struct  .u (nSteps x nInputs) physical input history
            %   Outputs:
            %       lb    (nParams x 1) double
            %       ub    (nParams x 1) double
            %       A     (nRows x nParams) double, sparse
            %       b     (nRows x 1) double
            arguments
                obj
                hist (1,1) struct
            end
            p     = obj.nInputs();
            L     = obj.nHorizon;
            uLast = obj.lastAppliedInput(hist, p);

            amp  = obj.constraints.amplitude;    % (p x 2)
            rate = obj.constraints.rate;         % (p x 2)

            % The trailing 1 is load-bearing: repelem of a SCALAR returns a row,
            % so the single-channel case would silently produce (1 x L) bounds.
            lb = repelem(amp(:, 1), L, 1);
            ub = repelem(amp(:, 2), L, 1);

            % The first sample of each channel additionally answers to u(k-1).
            first     = ((0:p-1) * L) + 1;
            lb(first) = max(lb(first), uLast(:) + rate(:, 1));
            ub(first) = min(ub(first), uLast(:) + rate(:, 2));

            % Never silently clamp: an empty box is reported, not repaired. It means
            % the last applied input lies further outside the amplitude box than one
            % rate step, which no admissible trajectory can recover from in one
            % sample -- a declaration or a history problem, not a rounding one.
            bad = find(lb(first) > ub(first), 1);
            if ~isempty(bad)
                error("rid:gen:emptyFeasibleSet", ...
                    "Channel ""%s"" has no feasible first sample: the actuator " + ...
                    "is at %g, outside the admissible amplitude range [%g, %g], " + ...
                    "and the rate limit [%g, %g] cannot bring it back inside in " + ...
                    "one sample -- the intersection is the empty interval " + ...
                    "[%g, %g]. Drive the actuator into the admissible range " + ...
                    "before starting the design, or widen the declaration.", ...
                    obj.constraints.channelName(bad), uLast(bad), ...
                    amp(bad, 1), amp(bad, 2), rate(bad, 1), rate(bad, 2), ...
                    lb(first(bad)), ub(first(bad)));
            end

            % Interior rate constraints. At L = 1 there are none and D is empty,
            % which is correct rather than a special case: the only rate limit that
            % applies is the one already folded into lb/ub above.
            D  = spdiags([-ones(L-1, 1), ones(L-1, 1)], 0:1, L-1, L);
            Dp = kron(speye(p), D);
            A  = [Dp; -Dp];
            b  = [repelem(rate(:, 2), L-1, 1); -repelem(rate(:, 1), L-1, 1)];
        end

        function [Theta0, info] = initialGuesses(obj, hist, thetaPrev, opts)
            %INITIALGUESSES  The multi-start set: warm start plus target sequences.
            %
            %   `+gen` CONSTRUCTS the starting points, because they must be feasible
            %   and expressed in the parametrization's own coordinates; `+opt`
            %   decides how many to run and what to do with the results. `nStarts`
            %   is therefore a REQUIRED argument with no default here --
            %   the count is declared at `opts.opt.nStarts`, where the decision
            %   lives, and this method merely answers. The dead-on-arrival
            %   diagnostic likewise belongs to `+opt`: it is read from the
            %   objective's `info.nShell`, which this class never sees.
            %
            %   ROW COMPOSITION:
            %
            %       row 1, from the second step on   the WARM START, shift(thetaPrev)
            %       the last row, always             an UNGUIDED target sequence
            %       every other row                  a GUIDED target sequence
            %
            %   On the first step there is no previous solution, so every row is a
            %   target sequence and the warm-start row does not exist.
            %   `info.warmStart` records which case applied. There is deliberately
            %   no constant-hold start.
            %
            %   THE TARGET SEQUENCE. Per row and per channel, repeat until the
            %   horizon is full: draw a target level; draw a speed fraction
            %   `g ~ U(0, 1]`; each sample take
            %
            %       step = min( g * largest legal step toward the target,
            %                   remaining distance )
            %
            %   clipped into the feasible increment interval; on landing EXACTLY on
            %   the target, draw a new target and a new `g`.
            %
            %   Far from the target this is a constant-rate ramp at `g * duMax`, and
            %   the last step lands exactly, so THERE IS NO ARRIVAL TOLERANCE -- the
            %   landing sample is assigned the target itself, not a value near it.
            %   `g` is drawn per target rather than per row, so one row visits
            %   several levels at several rates and spreads the `du` coordinate
            %   within a single start. `g` is drawn as `1 - rand`, because `rand`
            %   returns [0, 1) and a `g` of exactly zero would freeze the row.
            %
            %   FEASIBILITY. Before each draw the legal increments are
            %
            %       lo = max(duMin, uMin - u)      hi = min(duMax, uMax - u)
            %
            %   -- the rate limit intersected with what the amplitude box still
            %   allows. Never empty, because `lo <= 0 <= hi`: the rate box contains
            %   zero (rid.gen.Constraints enforces it) and `u` is already inside the
            %   amplitude box. Every sample is legal in both senses by construction,
            %   so no rejection, projection or repair is needed.
            %
            %   Inputs:
            %       hist               (1x1) struct  .u (nSteps x nInputs)
            %       thetaPrev          (nParams x 1) previous solution; EMPTY on
            %                                        the first step
            %       opts.nStarts       (1x1)         rows of Theta0; REQUIRED,
            %                                        no default
            %       opts.seed          (1x1)         nonnegative integer, or a stream
            %       opts.targetLevels  (:,1)         guided target levels;
            %                                        absent => every row unguided.
            %                                        SISO only: see the guard below
            %       opts.roiInputRange (nInputs x 2) the ROI's own u range; absent
            %                                        => the amplitude box, recorded
            %   Outputs:
            %       Theta0  (nStarts x nParams) double
            %       info    (1x1) struct  .seed .nStarts .warmStart .nSequences
            %                             .guided .targetRange .nTargets
            arguments
                obj
                hist         (1,1) struct
                thetaPrev    (:,1) double {mustBeReal} = zeros(0, 1)
                % NO DEFAULT. The number of local searches is a
                % search-strategy setting and lives at opts.opt.nStarts: +opt
                % decides the count, +gen answers. Giving it a default here would
                % put a second declaration of one number in the package that does
                % not own the decision -- the defect the coordinate-meaning rule exists to prevent.
                % Declared without a default value, on the epsReq0 precedent
                %, so an omitted argument is absent from opts and can be
                % refused by name rather than defaulted by nargin.
                opts.nStarts (1,1) double {mustBeInteger, mustBePositive}
                opts.seed = rid.defaultOptions().gen.seed
                opts.targetLevels  (:,:) double {mustBeReal} = zeros(0, 0)
                opts.roiInputRange (:,2) double {mustBeReal} = zeros(0, 2)
            end
            if ~isfield(opts, "nStarts")
                error("rid:gen:missingStartCount", ...
                    "nStarts is required and has no default here: how many " + ...
                    "local searches to run is the optimizer's decision, not " + ...
                    "the parametrization's, and it is declared at " + ...
                    "opts.opt.nStarts. +opt decides the count; +gen " + ...
                    "answers.");
            end
            p     = obj.nInputs();
            L     = obj.nHorizon;
            uLast = obj.lastAppliedInput(hist, p);
            amp   = obj.constraints.amplitude;
            rate  = obj.constraints.rate;

            % targetLevels ENTRIES ARE INPUT COMBINATIONS, so with several input
            % channels only a whole ROW has a defined resting output: the steady
            % state depends on every channel jointly. Drawing per column would mix
            % components of different rows into a combination with no relationship
            % to the output grid the surrogate-guided initialization built, and drawing per row is not well
            % defined either as the construction stands, because channels reach
            % their targets at different times and would have to redraw
            % independently. At nInputs = 1 the two readings coincide, which is why
            % nothing showed up before this guard.
            if ~isempty(opts.targetLevels) && p > 1
                error("rid:gen:targetLevelsNotDefinedForMiso", ...
                    "targetLevels was supplied for %d input channels, but its " + ...
                    "entries are input COMBINATIONS: with several channels the " + ...
                    "resting output depends on all of them jointly, so only a " + ...
                    "whole row has a defined output and drawing per channel " + ...
                    "would mix components of different rows. Multi-input runs " + ...
                    "initialize UNGUIDED instead -- omit targetLevels. The " + ...
                    "synchronized-segment design is not implemented.", p);
            end
            if ~isempty(opts.targetLevels) && size(opts.targetLevels, 2) ~= p
                error("rid:gen:invalidTargetLevels", ...
                    "targetLevels has %d column(s); the representation declares " + ...
                    "%d input channel(s).", size(opts.targetLevels, 2), p);
            end
            if ~isempty(opts.roiInputRange) && size(opts.roiInputRange, 1) ~= p
                error("rid:gen:invalidRoiInputRange", ...
                    "roiInputRange has %d row(s); the representation declares " + ...
                    "%d input channel(s).", size(opts.roiInputRange, 1), p);
            end

            % The unguided row draws over the ROI's own u range. Without it the
            % amplitude box stands in, and info says so rather than the substitution
            % being invisible.
            if isempty(opts.roiInputRange)
                roiRange    = amp;
                targetRange = "amplitudeBox";
            else
                roiRange    = opts.roiInputRange;
                targetRange = "roiInputRange";
            end

            % a target from EITHER source is clipped into the ROI u range and
            % the amplitude box. The ROI declares its u and y bounds independently,
            % so a level mapped through the static map can fall outside; clipping
            % plus a reachabilityReport note is the agreed handling, never silent
            % absorption.
            clipLo = max(roiRange(:, 1), amp(:, 1));
            clipHi = min(roiRange(:, 2), amp(:, 2));

            stream = rid.gen.Direct.buildStream(opts.seed);

            hasWarm = ~isempty(thetaPrev);
            Theta0  = zeros(opts.nStarts, obj.nParams());
            guided  = false(1, opts.nStarts);
            nTargets = zeros(opts.nStarts, p);

            firstSeq = 1;
            if hasWarm
                Theta0(1, :) = obj.shift(thetaPrev).';
                firstSeq     = 2;
            end
            % Guided rows use targetLevels; the LAST row is always unguided, and
            % without targetLevels every row is.
            for iStart = firstSeq:opts.nStarts
                isGuided = ~isempty(opts.targetLevels) && iStart < opts.nStarts;
                guided(iStart) = isGuided;
                uSeq = zeros(L, p);
                for cIn = 1:p
                    if isGuided
                        levels = opts.targetLevels(:, cIn);
                    else
                        levels = [];        % continuous draw over roiRange
                    end
                    [uSeq(:, cIn), nTargets(iStart, cIn)] = ...
                        rid.gen.Direct.targetSequence(uLast(cIn), L, ...
                            amp(cIn, :), rate(cIn, :), ...
                            [clipLo(cIn), clipHi(cIn)], ...
                            [roiRange(cIn, 1), roiRange(cIn, 2)], levels, stream);
                end
                Theta0(iStart, :) = uSeq(:).';
            end

            % The warm-start initialization says no repair is ever needed. This makes that a GUARANTEE
            % rather than a claim: it should never fire, and when it does the
            % construction is wrong, not the start. It caught the reflection
            % failure that the first construction had.
            [lb, ub, A, b] = obj.feasibleSet(hist);
            for iStart = 1:opts.nStarts
                theta = Theta0(iStart, :).';
                viol  = max([0; lb - theta; theta - ub]);
                if ~isempty(b)
                    viol = max(viol, max([0; A * theta - b]));
                end
                if viol > 1e-9
                    error("rid:gen:infeasibleStart", ...
                        "Start %d of %d violates the feasible set by %g. The " + ...
                        "construction is wrong: the warm-start initialization starts are feasible by " + ...
                        "construction and there is no repair path.", ...
                        iStart, opts.nStarts, viol);
                end
            end

            info = struct("seed", opts.seed, "nStarts", opts.nStarts, ...
                          "warmStart", string(missing), "nSequences", 0, ...
                          "guided", guided, "targetRange", targetRange, ...
                          "nTargets", nTargets);
            if hasWarm, info.warmStart = "shift"; else, info.warmStart = "none"; end
            info.nSequences = opts.nStarts - double(hasWarm);
        end

        function [pen, gradPen, info] = penalty(obj, Zphys, opts)
            %PENALTY  The soft-constraint term of the OBJECTIVE.
            %
            %   objective(theta) = J(theta) + penalty(Zphys(theta))
            %
            %   `pen` is the WEIGHTED total, so the caller adds it directly. It is
            %   never a term in `J`: `J` is a measurement of distribution quality
            %   and must stay comparable across runs, and `+crit` sees no constraint
            %   of any kind. The run record logs the two separately.
            %
            %   PHYSICAL COORDINATES. `Zphys` is in the representation's physical
            %   units, because `gMid`, `gHalf` and `s` come from the physical
            %   admissible declaration and `v` is made dimensionless by dividing by
            %   the physical admissible width. Note that the `Zphys` of
            %   the method description is the NORMALIZED one (its step 4 is
            %   part of the composition), so the two are not the same array.
            %
            %   Per soft constraint i on a simulated quantity g_i (section 4):
            %
            %       v_i(j) = ( |g_i(j) - gMid_i| - gHalf_i ) / s_i
            %       p_i    = mean_j  gamma(j) * max(0, v_i(j))^2
            %       w_i    = psiMax / violationScale^2
            %
            %   The violation is measured from the CENTRE of the declared box, not
            %   from zero: a limit is a box and is not in general symmetric about the
            %   origin, and writing `|g| - c` would silently assume `gMin = -gMax`.
            %   The aggregation is a MEAN, so an isolated violation is diluted by
            %   `1/nHorizon` -- the table in section 4 reads that consequence out.
            %
            %   The weight is expressed through what it COSTS rather than as a raw
            %   number: `psiMax` is read from `opts.crit` and never duplicated under
            %   `opts.gen`, so the calibration cannot drift when `psiMax` changes
            %   (the psiMax = 1 default is open).
            %
            %   Inputs:
            %       Zphys                 (nHorizon x nDim) double, PHYSICAL
            %       opts.violationScale  (1x1) fraction of the admissible width whose
            %                            sustained violation costs psiMax
            %       opts.tailWeight      (1x1) gamma(end)/gamma(1)
            %       opts.psiMax          (1x1) from opts.crit, never from opts.gen
            %   Outputs:
            %       pen      (1x1) double                weighted total
            %       gradPen  (nHorizon x nDim) double    d(pen)/d(Zphys), NOT d/dtheta
            %       info     (1x1) struct                .w .gamma .p .name .nActive
            %
            %   `gradPen` is with respect to `Zphys` and not to `theta` on purpose: it
            %   then shares the simulation Jacobian with the criterion's gradient
            %   path instead of recomputing it.
            arguments
                obj
                Zphys (:,:) double {mustBeReal}
                opts.violationScale (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().gen.violationScale
                opts.tailWeight     (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().gen.tailWeight
                opts.psiMax         (1,1) double {mustBePositive} ...
                    = rid.defaultOptions().crit.psiMax
            end
            L    = obj.nHorizon;
            nDim = obj.constraints.nDim;
            if size(Zphys, 1) ~= L
                error("rid:gen:horizonMismatch", ...
                    "Zphys has %d row(s); the horizon is %d.", size(Zphys, 1), L);
            end
            if size(Zphys, 2) ~= nDim
                error("rid:gen:invalidPrediction", ...
                    "Zphys has %d column(s); the representation has %d " + ...
                    "dimension(s).", size(Zphys, 2), nDim);
            end

            gamma = rid.gen.Direct.horizonWeights(L, opts.tailWeight);
            w     = opts.psiMax / opts.violationScale^2;
            soft  = obj.constraints.soft;

            pen     = 0;
            gradPen = zeros(L, nDim);
            pPer    = zeros(1, numel(soft));
            for iSoft = 1:numel(soft)
                thisSoft = soft(iSoft);
                offset   = Zphys(:, thisSoft.iDim) - thisSoft.gMid;
                v        = (abs(offset) - thisSoft.gHalf) / thisSoft.s;
                vPos     = max(0, v);

                pPer(iSoft) = mean(gamma .* vPos.^2);
                pen         = pen + w * pPer(iSoft);

                % d/dg of mean_j gamma_j max(0,v_j)^2, with dv/dg = sign(g-gMid)/s.
                % At g == gMid the sign is 0, but there vPos is 0 as well, so the
                % term vanishes either way.
                gradPen(:, thisSoft.iDim) = gradPen(:, thisSoft.iDim) + ...
                    w * gamma .* (2 * vPos) .* sign(offset) / (thisSoft.s * L);
            end

            info = struct("w", w, "gamma", gamma, "p", pPer, ...
                          "name", [soft.name], "nActive", sum(pPer > 0));
        end

        function s = describe(obj)
            %DESCRIBE  Human-readable summary, for logs and run records.
            %
            %   Outputs:
            %       s  (1x1) string
            s = sprintf("rid.gen.Direct  nHorizon = %d, nInputs = %d, " + ...
                "nParams = %d\n  applied unit : 1 sample\n" + ...
                "  theta layout : channel-major, uSeq(:)\n%s", ...
                obj.nHorizon, obj.nInputs(), obj.nParams(), ...
                obj.constraints.describe());
        end

    end

    methods (Static)

        function [uSeq, nTargets] = targetSequence(u0, L, amp, rate, clipBox, ...
                                                    roiRange, levels, stream)
            %TARGETSEQUENCE  One channel of one start: ramp between drawn targets.
            %
            %   `levels` empty means draw the target continuously over `roiRange`
            %   (the unguided row); otherwise draw uniformly from `levels`.
            %   Either way the target is clipped into `clipBox`, the intersection of
            %   the ROI u range with the amplitude box.
            %
            %   Inputs:
            %       u0        (1x1)  u(k-1), the level the channel starts from
            %       L         (1x1)  horizon
            %       amp       (1x2)  amplitude box of this channel
            %       rate      (1x2)  rate limit of this channel
            %       clipBox   (1x2)  where a target is allowed to be
            %       roiRange  (1x2)  where an unguided target is drawn from
            %       levels    (:,1)  guided levels, or empty
            %       stream           RandStream, or any object answering rand
            %   Outputs:
            %       uSeq      (L x 1) double
            %       nTargets  (1x1)   how many targets the row visited
            uSeq     = zeros(L, 1);
            u        = u0;
            jStep    = 0;
            nTargets = 0;

            while jStep < L
                if isempty(levels)
                    target = roiRange(1) + (roiRange(2) - roiRange(1)) * rand(stream);
                else
                    target = levels(min(numel(levels), ...
                                        1 + floor(numel(levels) * rand(stream))));
                end
                target   = min(max(target, clipBox(1)), clipBox(2));
                g        = 1 - rand(stream);      % rand is [0,1); g must be (0,1]
                nTargets = nTargets + 1;

                if target == u
                    % Already there. Consume one sample holding still -- legal,
                    % because the rate box contains zero -- so the loop cannot spin
                    % forever on a degenerate level set.
                    jStep        = jStep + 1;
                    uSeq(jStep)  = u;
                    continue
                end

                while jStep < L
                    lo   = max(rate(1), amp(1) - u);
                    hi   = min(rate(2), amp(2) - u);
                    gap  = target - u;
                    if gap > 0, largest = hi; else, largest = -lo; end
                    stepMag = min(g * largest, abs(gap));
                    jStep   = jStep + 1;
                    if stepMag >= abs(gap)
                        % The landing sample IS the target, assigned rather than
                        % accumulated: there is no arrival tolerance to choose, and
                        % bit equality is what the section 10 test asserts.
                        u = target;
                    else
                        step = min(max(sign(gap) * stepMag, lo), hi);
                        u    = u + step;
                    end
                    uSeq(jStep) = u;
                    if u == target, break, end
                end
            end
        end

        function stream = buildStream(seed)
            %BUILDSTREAM  A private stream from a seed or an existing RandStream.
            %
            %   The global stream is never read or written, so a start set
            %   is reproducible in isolation and inside parfor. A non-numeric seed
            %   is passed through and used duck-typed via `rand(stream)`, so a test
            %   can script the draws without the production path gaining a knob.
            if ~isnumeric(seed)
                stream = seed;
            else
                stream = RandStream("twister", "Seed", seed);
            end
        end

        function gamma = horizonWeights(nHorizon, tailWeight)
            %HORIZONWEIGHTS  The horizon weighting gamma, normalized to UNIT MEAN.
            %
            %       gammaRaw(j) = tailWeight ^ ( (j-1) / (nHorizon-1) )
            %       gamma(j)    = gammaRaw(j) / mean_j gammaRaw(j)
            %
            %   Geometric IN THE EXPONENT, so the profile is horizon-length
            %   invariant: the far end always carries `tailWeight` of the applied
            %   sample's weight whatever `nHorizon` is, where a fixed decay rate
            %   would silently change the shape when the horizon changed.
            %
            %   THE UNIT-MEAN NORMALIZATION IS WHAT MAKES gamma A SHAPE RATHER THAN
            %   ALSO A MAGNITUDE. Without it `mean(gammaRaw) < 1` -- 0.399 at
            %   tailWeight = 0.1, nHorizon = 20 -- and the whole `violationScale`
            %   calibration would be wrong by that factor. Normalized, `tailWeight`
            %   decides only WHERE in the horizon a violation hurts and
            %   `violationScale` decides HOW MUCH, so the two knobs cannot interfere.
            %
            %   Inputs:
            %       nHorizon    (1x1) double  positive integer
            %       tailWeight  (1x1) double  positive
            %   Outputs:
            %       gamma       (nHorizon x 1) double, mean exactly 1
            arguments
                nHorizon   (1,1) double {mustBeInteger, mustBePositive}
                tailWeight (1,1) double {mustBePositive}
            end
            if nHorizon == 1
                % A single position carries the whole horizon: there is no ratio to
                % express, and (j-1)/(nHorizon-1) would be 0/0.
                gamma = 1;
                return
            end
            raw   = tailWeight .^ (((1:nHorizon).' - 1) / (nHorizon - 1));
            gamma = raw / mean(raw);
        end

    end

    methods (Access = private)

        function assertTrajectorySize(obj, uSeq)
            %ASSERTTRAJECTORYSIZE  A trajectory is (nHorizon x nInputs), always.
            if ~isequal(size(uSeq), [obj.nHorizon, obj.nInputs()])
                error("rid:gen:invalidTrajectory", ...
                    "uSeq is %s; this parametrization produces (%d x %d).", ...
                    mat2str(size(uSeq)), obj.nHorizon, obj.nInputs());
            end
        end

        function uLast = lastAppliedInput(obj, hist, p)
            %LASTAPPLIEDINPUT  u(k-1), the state the feasible set depends on.
            %
            %   REQUIRED, NEVER DEFAULTED, on the precedent rather than the
            %   warm-up row rule one. The warm-up row rule is about warm-up rows, which the start of a run
            %   genuinely does not have; this is not that. A real actuator sits at
            %   some physical position when the experiment begins -- the value is
            %   UNRECORDED, not undefined -- and the xStart requirement already requires the run to
            %   know it, since steadyState(uConst) needs it to produce xStart.
            %
            %   Defaulting it would leave the very FIRST commanded step unbounded by
            %   any rate limit: the step most likely to damage equipment, and the one
            %   no later step can undo.
            % Structural first; what follows is the SEMANTIC question,
            % and its message is the reason it stays here: it tells a user what to
            % do, which a generic identifier would throw away.
            rid.assertHistory(hist);
            if isempty(hist.u)
                error("rid:gen:missingInitialInput", ...
                    "feasibleSet needs the input the process starts from, as the " + ...
                    "last row of hist.u; it is absent or empty. Declare the input " + ...
                    "the process starts from -- it is unrecorded, not undefined, " + ...
                    "and without it the first commanded step carries no rate limit.");
            end
            if size(hist.u, 2) ~= p
                error("rid:gen:invalidHistory", ...
                    "hist.u has %d column(s); the representation declares %d " + ...
                    "input channel(s).", size(hist.u, 2), p);
            end
            uLast = hist.u(end, :);
        end

    end

end
