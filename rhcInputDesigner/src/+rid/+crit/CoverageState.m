classdef CoverageState < rid.crit.Criterion
%COVERAGESTATE  Incremental evaluation of the criterion, with a raw-distance cache.
%
%   THE FIRST IMPLEMENTATION OF rid.crit.Criterion, and the one that
%   overrides every default there: it certifies (`isSatisfied`), is monotone under
%   `add` (invariant 3), has a demanded resolution, an active set and supporting
%   points. Its behaviour did not change when the interface was introduced; the
%   four stored figure runs replay bit-exactly in u and in the J trace.
%
%   The criterion sits in the innermost loop of the receding-horizon optimization,
%   so scoring a candidate must not cost a pass over the whole history. This class
%   holds the nearest-neighbour distance of every supporting point to the committed
%   samples and updates it in place:
%
%       add(Znew)      dRaw = min(dRaw, distances to Znew)   O(nSup * nNew)
%       probe(Zcand)   score without committing               O(nActive * nCand)
%       value()        the current J                          O(nSup)
%       isSatisfied()  has the criterion certified itself?    O(nSup)
%
%   THE CACHE HOLDS THE RAW, UNNORMALIZED DISTANCE -- never dTilde. For a fixed
%   supporting point, epsReq(j) is a constant scale factor, so argmin_i is
%   unaffected by the normalization: caching dRaw is exactly equivalent and stays
%   valid when epsReq changes, which it does whenever the refinement fields are
%   re-estimated. Anything cached in units of epsReq goes stale at that moment,
%   silently. `setResolution` is therefore O(nSup) and does not touch the cache at
%   all (section 7.1).
%
%   THE ACTIVE SET: ADD MAINTAINS IT, PROBE EXPLOITS IT. The class tracks
%   which supporting points are not yet resolved, and a candidate evaluation
%   measures distances only to those -- O(nActive * nHorizon).
%
%   Why skipping is exact in `probe`: if `dRaw_j <= epsReq_j` then
%   `min(dRaw_j, anything) <= epsReq_j` too, so `psi_j = 0` and `isShell_j` is
%   false for EVERY candidate. A resolved point cannot be made to contribute by a
%   candidate, to `J` or to the gradient, so leaving it out changes no number.
%
%   Why `add` measures EVERYTHING anyway, though a resolved point cannot come back
%   on its own: because `dRaw` must stay TRUE, not merely a bound. Let `epsReq = 1`
%   and take three operations -- `add` at distance 0.9 (resolved), `add` at
%   distance 0.5, `setResolution(0.6)`. With a skipping `add`, `dRaw` is frozen at
%   0.9, the point un-resolves at `dTilde = 1.5` and `J = 0.5`; the truth is
%   `dRaw = 0.5`, `dTilde = 0.833` and `J = 0`. The frozen value is an upper bound,
%   which is invisible while `epsReq` holds and wrong the moment it moves. `add`
%   costs `O(nSup * nNew)` with `nNew = 1` per loop step -- nothing, against
%   `nStarts * iterations` probes -- so the saving was never there to take.
%
%   SATURATED POINTS DO NOT LEAVE, and this is the tempting mistake. A point beyond
%   `1 + psiMax` also contributes zero gradient and a constant to `J`, so it looks
%   equally skippable. It is not: `dRaw` only decreases, so a saturated point can
%   move INTO the shell and start contributing. The criterion for leaving is
%   MONOTONICITY, not zero gradient.
%
%   TWO INVALIDATIONS, both on paths that already existed. `setResolution`: a
%   rising learned complexity shrinks `epsReq_j`, raises `dTilde_j` and can
%   un-resolve a point -- this is invariant 2, and "resolved" is exactly a quantity
%   held in units of `epsReq`. `rebuild`: a surrogate update moves historical
%   coordinates (invariant 3) and `dRaw` can come back larger. Membership is
%   recomputed at both, riding on the invalidation the cache already has.
%
%   BOOKKEEPING: `J` is a mean over ALL nSup points. Resolved points contribute
%   zero to the numerator and one to the denominator. The active set skips
%   COMPUTATION, not membership; averaging over it alone would inflate `J`.
%
%   In certification mode `J = 0` means every point is resolved, so the set empties
%   exactly at termination: the design gets cheaper the closer it is to finishing,
%   and the expensive part of a run is the part with most left to gain.
%
%   The complexity contract -- O(nActive * nHorizon) per candidate, never O(nSup * N)
%   -- IS TESTED, by counting. Every distance this class computes passes through
%   one protected seam, `distanceMatrix`, which a test subclass overrides to count
%   `numel`; CoverageStateGradientTest asserts the exact count against a history
%   an order of magnitude longer than the candidate block, and
%   CoverageStateActiveSetTest asserts that the count per probe FALLS as points
%   resolve while the count per `add` stays flat -- which is what pins the split
%   between the two.
%
%   This reverses what this header said earlier, and the reversal has two
%   reasons. The loop's inner iteration now depends on algorithmic invariant 1,
%   where a violation is INVISIBLE -- the run is merely slow, not wrong, and slow
%   for no reason a plot would show. And review-only enforcement has already failed
%   in this project once, on the golden layer, which was decorative on every
%   machine but one for its whole existence before anyone noticed.
%
%   The objection that a counting double changes the code to suit the test is
%   weaker than it looks here: collapsing two distance call sites into one named
%   seam is a simplification on its own terms, and it is the place the anisotropic
%   extension would have to change anyway.
%
%   The behavioural half is still asserted and still necessary: `probe` returns the
%   same value whether the cached dRaw came from ten committed samples or ten
%   thousand, because it reads only dRaw and never the history. Neither test
%   subsumes the other -- an implementation that computed all N distances and threw
%   them away would pass the behavioural test while violating the invariant, and
%   only the count catches it.
%
%   Correctness is asserted differently. The cached path has one failure mode -- a
%   stale cache -- and it is silent: no error, only slightly worse decisions,
%   indistinguishable from the criterion simply performing badly. So this class is
%   compared against rid.crit.evaluate over randomized sequences of
%   add / probe / setResolution / rebuild, to 1e-12. For that oracle to
%   mean anything the two must not share code, so the aggregation below is written
%   out again rather than factored into a helper both call.
%
%   Value semantics do not apply: this is the one place in the design, besides the
%   update schedulers, where mutable state is intended.
%   `probe` is nevertheless side-effect free, and the tests assert that on the
%   handle itself.
%
%   `probeGradient` is the cached counterpart of `rid.crit.gradient`, required by the optimizer: the reference path takes its value
%   from `rid.crit.evaluate`, which walks every sample, and nSamples grows for the
%   whole run. Wiring that into the optimizer would violate invariant 1.
%
%   Construction:
%       state = rid.crit.CoverageState(Q, epsReq)
%       state = rid.crit.CoverageState(Q, epsReq, "psiMax", 1, "pNorm", 1)
%
%   Inputs:
%       Q             (nSup x nDim)  supporting points, normalized coordinates
%       epsReq        (nSup x 1)     demanded resolution, finite and > 0
%       opts.psiMax   (1x1)          > 0; Inf allowed
%       opts.pNorm    (1x1)          >= 1; Inf allowed

    properties (SetAccess = immutable)
        nSup   (1,1) double   % the TRUE count, kept even when Q is removed
        nDim   (1,1) double
        psiMax (1,1) double   % fixed for the run: it sets the criterion's character
        pNorm  (1,1) double   % fixed for the run
    end

    properties (SetAccess = private)
        % (nSup x nDim) supporting points, FIXED FOR THE RUN. Private rather than
        % immutable for exactly one reason, and it is written down: `withoutQ`
        % has to be able to clear it, because `opts.io.storeQ = false` must
        % remove `Q` from a stored record completely and the class that owns the
        % invariant is the only thing that can produce a Q-less instance (the storeQ rule,
        % the single-setup rule). Nothing else writes it after construction.
        Q      (:,:) double
    end

    properties (SetAccess = private)
        epsReq   (:,1) double  % (nSup x 1) changed only by setResolution
        dRaw     (:,1) double  % (nSup x 1) THE CACHE: raw NN distance, Inf where uncovered
        isActive (:,1) logical % (nSup x 1) NOT yet resolved; probe measures only these
        nSamples (1,1) double  % committed samples behind dRaw; diagnostics only
    end

    methods

        function obj = CoverageState(Q, epsReq, opts)
            arguments
                Q             (:,:) double {mustBeReal, mustBeFinite}
                epsReq        (:,1) double {mustBeReal}
                opts.psiMax   (1,1) double {mustBeReal} = rid.defaultOptions().crit.psiMax
                opts.pNorm    (1,1) double {mustBeReal} = rid.defaultOptions().crit.pNorm
            end

            if size(Q, 1) < 1
                error("rid:crit:invalidResolution", ...
                    "Q must have at least one supporting point.");
            end
            if isnan(opts.psiMax) || opts.psiMax <= 0
                error("rid:crit:invalidPenaltyClip", ...
                    "psiMax must be strictly positive (Inf is allowed); got %g.", opts.psiMax);
            end
            if isnan(opts.pNorm) || opts.pNorm < 1
                error("rid:crit:invalidNorm", ...
                    "pNorm must be >= 1 (Inf is allowed); got %g.", opts.pNorm);
            end

            obj.Q      = Q;
            obj.nSup   = size(Q, 1);
            obj.nDim   = size(Q, 2);
            obj.psiMax = opts.psiMax;
            obj.pNorm  = opts.pNorm;

            obj.assertResolution(epsReq);
            obj.epsReq = epsReq;

            % Nothing committed yet: every supporting point is infinitely far from
            % a distribution with no samples, hence unresolved and active.
            obj.dRaw     = inf(obj.nSup, 1);
            obj.nSamples = 0;
            obj.refreshActive();
        end

        function add(obj, Znew)
            %ADD  Commit accepted samples and update the cache.
            %
            %   The nearest neighbour of a supporting point can only get nearer as
            %   samples are added, so the update is a running minimum -- exact, and
            %   independent of the order in which samples arrive.
            %
            %   ADD MAINTAINS THE ACTIVE SET; it does not exploit it. Every
            %   supporting point is measured, resolved or not, because `dRaw` must
            %   stay the TRUE nearest distance rather than an upper bound on it --
            %   see the three-operation counter-example in the class header. The
            %   saving would have been O(nSup) per step, which is nothing.
            %
            %   Inputs:
            %       Znew  (nNew x nDim)  committed samples; empty is a no-op
            arguments
                obj
                Znew (:,:) double {mustBeReal, mustBeFinite}
            end
            if isempty(Znew)
                return
            end
            obj.assertColumns(Znew, "Znew");
            obj.dRaw     = min(obj.dRaw, obj.nearestDistance(Znew));
            obj.nSamples = obj.nSamples + size(Znew, 1);
            obj.refreshActive();
        end

        function J = probe(obj, Zcand)
            %PROBE  Score a candidate without committing it.
            %
            %   Reads the cache and never writes it: this is the call the optimizer
            %   makes once per candidate, and a mutation here would make the
            %   candidates depend on the order they were scored in.
            %
            %   PROBE EXPLOITS THE ACTIVE SET: distances are measured only
            %   to unresolved supporting points, so this costs O(nActive * nCand).
            %   That is exact rather than approximate -- a resolved point has
            %   `dRaw <= epsReq`, so `min(dRaw, anything) <= epsReq` and `psi = 0`
            %   whatever the candidate is. The aggregation below still runs over all
            %   nSup rows: this skips computation, not membership.
            %
            %   Inputs:
            %       Zcand  (nCand x nDim)  candidate samples; empty scores the
            %                              current state
            %   Outputs:
            %       J      (1x1)
            arguments
                obj
                Zcand (:,:) double {mustBeReal, mustBeFinite}
            end
            if ~isempty(Zcand)
                obj.assertColumns(Zcand, "Zcand");
            end
            J = obj.aggregate(obj.candidateDistance(Zcand));
        end

        function [gradCand, info] = probeGradient(obj, Zcand)
            %PROBEGRADIENT  dJ/dZcand at the probed state, from the cache alone.
            %
            %   Returns one row per CANDIDATE sample -- never per history sample --
            %   so it lines up directly with the `A` of the affine map, which
            %   carries the horizon rows only.
            %   This class still sees nothing but points: the composition is the
            %   caller's, and naming the package that performs it here would be a
            %   dependency +crit does not have.
            %
            %   WHY THIS COSTS O(nSup * nHorizon) AND NOT O(nSup * N). A supporting
            %   point whose nearest sample lies in the realized history contributes
            %   EXACTLY ZERO to the gradient: that sample does not move with theta,
            %   so its Jacobian block is zero however deep inside the shell the
            %   point sits (section 3). The gradient therefore never needs a history
            %   sample's coordinates at all -- only the cached raw distance, to
            %   decide whether a candidate beats it. Per supporting point: compare
            %   the cached distance against the nHorizon candidate distances; a term
            %   exists only where a candidate WINS and the point lands in the shell.
            %
            %   That history result is a REQUIREMENT here, not a corollary, because
            %   a finite-difference reference makes the same restriction
            %   automatically: an implementation that got it wrong would agree with
            %   its own check, and only a constructed case exposes it (section 3).
            %
            %   TIES GO TO THE HISTORY, by the strict `<` below. That is not a free
            %   choice: `rid.crit.gradient` breaks ties by lowest row index, and the
            %   rows arrive in time order -- history first, horizon last
            %   -- so a history sample tied
            %   with a candidate wins there and contributes nothing.
            %   Ties among candidates go to the lowest candidate index, which `min`
            %   already does and which is the same row in either ordering. Both are
            %   genuine elements of the Clarke generalized gradient (section 5).
            %
            %   NO MUTATION, like `probe`: this is a scoring call, and only `add`
            %   commits.
            %
            %   The aggregation and the accumulation are written out again rather
            %   than shared with rid.crit.gradient, exactly as `aggregate` is not
            %   shared with rid.crit.evaluate: an oracle on the same code path
            %   checks nothing. The equivalence is asserted instead.
            %
            %   Inputs:
            %       Zcand     (nCand x nDim)  candidate samples; empty is legal and
            %                                 yields a zero-row gradient
            %   Outputs:
            %       gradCand  (nCand x nDim)  dJ/d(Zcand(i,:))
            %       info      (1x1) struct  .J .nearest (nSup x 1, the winning
            %                               CANDIDATE row, 0 where the history wins)
            %                               .fromCandidate (nSup x 1) .isShell
            %                               .nShell .dRaw .dTilde .psi
            arguments
                obj
                Zcand (:,:) double {mustBeReal, mustBeFinite}
            end

            % The p -> Inf limit is a MAX, so its derivative is that of the argmax
            % term alone -- a different expression, not this one evaluated at
            % pNorm = Inf. See rid.crit.gradient, which refuses it for the same
            % reason and states what the correct form is.
            if isinf(obj.pNorm)
                error("rid:crit:gradientNormUnsupported", ...
                    "The gradient is not implemented for pNorm = Inf. At the " + ...
                    "worst-case limit J is a max over supporting points, so the " + ...
                    "derivative is that of the argmax term alone; evaluating the " + ...
                    "p-norm expression there gives 1/nSup where 1 is correct, " + ...
                    "which is wrong by a factor of nSup and looks plausible.");
            end

            nCand    = size(Zcand, 1);
            gradCand = zeros(nCand, obj.nDim);

            nearest  = zeros(obj.nSup, 1);
            fromCand = false(obj.nSup, 1);
            dRawCand = obj.dRaw;

            % LIKE probe, THIS EXPLOITS THE ACTIVE SET and for the same
            % reason: a resolved point has dTilde <= 1 for every candidate, so it
            % is never in the shell and never contributes a gradient term.
            act = obj.isActive;
            if nCand > 0 && any(act)
                obj.assertColumns(Zcand, "Zcand");
                [dCand, iCand] = min(obj.distanceMatrix(Zcand, act), [], 2);
                % STRICT: a candidate must BEAT the cached distance, not merely
                % match it. See the tie paragraph in the header.
                wins            = false(obj.nSup, 1);
                wins(act)       = dCand < obj.dRaw(act);
                dRawCand(wins)  = dCand(wins(act));
                nearest(wins)   = iCand(wins(act));
                fromCand        = wins;
            elseif nCand > 0
                obj.assertColumns(Zcand, "Zcand");
            end

            dTilde = dRawCand ./ obj.epsReq;
            psi    = min(max(dTilde - 1, 0), obj.psiMax);
            J      = mean(psi.^obj.pNorm)^(1 / obj.pNorm);

            % Strict inequalities: a point exactly at a clip boundary is outside the
            % shell and contributes zero -- one of the two one-sided values, taken
            % consistently, and the same convention rid.crit.gradient uses.
            isShell = dTilde > 1 & dTilde < 1 + obj.psiMax;

            info = struct("J", J, "nearest", nearest, "fromCandidate", fromCand, ...
                          "isShell", isShell, "nShell", sum(isShell), ...
                          "dRaw", dRawCand, "dTilde", dTilde, "psi", psi);

            % J == 0 on EXACT equality, no tolerance: J is identically zero on a
            % neighbourhood unless some dTilde sits exactly at 1, so zero is the
            % TRUE gradient here and not a substitute (section 5).
            if J == 0 || nCand == 0
                return
            end

            for jSup = 1:obj.nSup
                if ~isShell(jSup) || ~fromCand(jSup)
                    continue
                end
                dRawJ = dRawCand(jSup);
                if dRawJ == 0
                    % uvec is 0/0; psi is 0 here anyway, so the term is doubly inert.
                    continue
                end
                iNear  = nearest(jSup);
                uvec   = (Zcand(iNear, :) - obj.Q(jSup, :)) / dRawJ;
                dJdPsi = (1 / obj.nSup) * (psi(jSup) / J)^(obj.pNorm - 1);
                gradCand(iNear, :) = gradCand(iNear, :) + ...
                    dJdPsi * (1 / obj.epsReq(jSup)) * uvec;
            end
        end

        function J = value(obj)
            %VALUE  The criterion at the committed samples.
            %
            %   Outputs:
            %       J  (1x1)
            J = obj.aggregate(obj.dRaw);
        end

        function tf = isSatisfied(obj)
            %ISSATISFIED  Has this criterion certified its own completion?
            %
            %   `J == 0`: every supporting point is within its demanded resolution,
            %   which is the certification condition of the method description
            %   and the termination condition of the method description.
            %
            %   IT EXISTS SO THE LOOP DOES NOT TEST `J == 0` ITSELF. That
            %   test is a statement about what certification MEANS, and
            %   the method description keeps criterion maths out of
            %   `+rhc`. It is `false` by default on any other criterion: `J = 0`
            %   certification is a property of THIS one, not of criteria in general.
            %
            %   READ FROM `value()`, NOT FROM THE ACTIVE SET. `nActive == 0` and
            %   `J == 0` are equivalent, and section 6 says the run stops on `J` so
            %   that one fact has one source -- `nActive` is logged and never
            %   terminates. Routing this through `value()` keeps that true.
            %
            %   EXACT EQUALITY, no tolerance, for the reason
            %   the method description gives for the gradient's zero
            %   branch: `max(dTilde - 1, 0)` returns an exact zero, so `J` is
            %   exactly zero when every point is resolved, and a `J < tol` test
            %   would certify a run that has not finished.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = obj.value() == 0;
        end

        function tf = isMonotone(obj) %#ok<MANU>
            %ISMONOTONE  Yes: `add` only shrinks dRaw, so J is non-increasing (P2).
            %
            %   Holds between `rebuild` and `setResolution` only; the loop resets
            %   its monotonicity baseline at both (invariant 3).
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = true;
        end

        function epsReq = resolution(obj)
            %RESOLUTION  The demanded resolution in force, (nSup x 1).
            %
            %   Outputs:
            %       epsReq  (nSup x 1) double
            epsReq = obj.epsReq;
        end

        function n = nActive(obj)
            %NACTIVE  Supporting points not yet resolved.
            %
            %   Outputs:
            %       n  (1x1) double
            n = sum(obj.isActive);
        end

        function tf = isConsistentWith(obj, Q)
            %ISCONSISTENTWITH  Was this state constructed against exactly this Q?
            %
            %   The cache holds distances to the state's OWN supporting points, so
            %   a different Q would have the loop report coverage of one set while
            %   the record describes another (rid.rhc.Setup).
            %
            %   Inputs:
            %       Q   (nSup x nDim) double
            %   Outputs:
            %       tf  (1x1) logical
            tf = isequal(obj.Q, Q);
        end

        function fresh = freshCopy(obj)
            %FRESHCOPY  Same Q, resolution in force, psiMax and pNorm; nothing committed.
            %
            %   NOT what `rid.io.replay` uses for this class: a replay must start
            %   from the resolution the run STARTED with, which the record holds
            %   as `epsReqInitial`, and this object holds the one it
            %   ended with.
            %
            %   Outputs:
            %       fresh  (1x1) rid.crit.CoverageState
            fresh = rid.crit.CoverageState(obj.Q, obj.epsReq, ...
                        "psiMax", obj.psiMax, "pNorm", obj.pNorm);
        end

        function s = describe(obj)
            %DESCRIBE  One line for the run record.
            %
            %   Outputs:
            %       s  (1x1) string
            s = sprintf("rid.crit.CoverageState: %d supporting point(s) in %d " + ...
                "dimension(s), psiMax %g, pNorm %g", obj.nSup, obj.nDim, ...
                obj.psiMax, obj.pNorm);
        end

        function setResolution(obj, epsReq)
            %SETRESOLUTION  Adopt a new demanded resolution, in O(nSup).
            %
            %   The cache is deliberately untouched: it holds raw distances, which
            %   do not depend on epsReq. This is the
            %   call made after a refinement-field update.
            %
            %   THE ACTIVE SET IS RECOMPUTED HERE, and this is invariant 2 acting:
            %   "resolved" is a quantity held in units of epsReq, so a rising
            %   learned complexity that shrinks epsReq_j raises dTilde_j and can
            %   un-resolve a point. The cache itself is untouched and stays exact,
            %   which is what makes the un-resolved point's dTilde correct rather
            %   than an overestimate.
            %
            %   Inputs:
            %       epsReq  (nSup x 1)  finite and strictly positive
            arguments
                obj
                epsReq (:,1) double {mustBeReal}
            end
            obj.assertResolution(epsReq);
            obj.epsReq = epsReq;
            obj.refreshActive();
        end

        function rebuild(obj, Zall)
            %REBUILD  Recompute the cache from the full history.
            %
            %   Needed when past coordinates MOVE -- with simulated past values, or
            %   any representation with dependsOnSurrogate = true, a surrogate
            %   update relocates every point in regressor space and the incremental
            %   minimum no longer describes the same geometry (section 7.3). The
            %   caller reads that condition from rep.dependsOnSurrogate and
            %   pastSource; this class is told, never assumes.
            %
            %   AN EMPTY `Zall` IS THE RESET, and it is used deliberately: the
            %   budget-mode coverage reset of the coverage reset clears the cache and returns
            %   every supporting point to the active set, so the design proceeds as
            %   though no data existed. That is a fiction the caller chooses and
            %   logs; this class merely performs it, and the honest coverage is
            %   recovered by one `rebuild(Zall)` over the real history at run end.
            %
            %   Inputs:
            %       Zall  (nAll x nDim)  the full history in current coordinates;
            %                            empty resets to the uncovered state
            arguments
                obj
                Zall (:,:) double {mustBeReal, mustBeFinite}
            end
            if isempty(Zall)
                obj.dRaw     = inf(obj.nSup, 1);
                obj.nSamples = 0;
                obj.refreshActive();
                return
            end
            obj.assertColumns(Zall, "Zall");
            obj.dRaw     = obj.nearestDistance(Zall);
            obj.nSamples = size(Zall, 1);
            % A surrogate update moves historical coordinates (invariant 3), so
            % dRaw can come back LARGER and a resolved point can un-resolve.
            obj.refreshActive();
        end

        function copy = withoutQ(obj)
            %WITHOUTQ  A COPY of this state with the supporting points removed.
            %
            %   FOR A RUN RECORD, and for nothing else. `storeQ = false`
            %   must remove `Q` from the record completely -- `problem.Q`,
            %   `setup.Q` and this object alike -- because the flag exists for
            %   SIZE, and while `Q` sat inside the stored `Setup` the file never
            %   actually got smaller. A record could then announce
            %   `storedQ = false` and carry the points anyway: a record
            %   misreporting itself.
            %
            %   IT RETURNS A COPY, and that is not a detail. This is the
            %   project's one intentional handle class, so clearing `Q` in place
            %   would empty the live state the caller is still holding -- the
            %   run's own `result` losing its supporting points as a side effect
            %   of being recorded.
            %
            %   EVERYTHING ELSE SURVIVES: `epsReq`, the distance cache, the
            %   active set, `nSup`, `nDim`, `psiMax` and `pNorm`. So
            %   `rid.crit.coverageReport` still answers on a stripped record,
            %   which is the whole point of removing only the coordinates. In
            %   particular `nSup` stays the TRUE count: a record reporting
            %   `nSup = 0` would describe a run that had no supporting points at
            %   all, a different and worse untruth than one that omits where
            %   they were.
            %
            %   Outputs:
            %       copy  (1x1) rid.crit.CoverageState  with `Q` (0 x nDim)
            copy = rid.crit.CoverageState(obj.Q, obj.epsReq, ...
                       "psiMax", obj.psiMax, "pNorm", obj.pNorm);

            % Carried over, not recomputed: a rebuilt cache would silently reset
            % dRaw to Inf and report a run that had covered nothing.
            copy.dRaw     = obj.dRaw;
            copy.isActive = obj.isActive;
            copy.nSamples = obj.nSamples;

            copy.Q = zeros(0, obj.nDim);
        end

    end

    methods (Access = protected)

        function D = distanceMatrix(obj, Znorm, rows)
            %DISTANCEMATRIX  Every distance this class computes, (nnz(rows) x nRows).
            %
            %   THE SINGLE DISTANCE SEAM. `nearestDistance` needs the minimum and
            %   `probeGradient` needs the argmin as well, and routing both through
            %   one method means the class has exactly one place where distance
            %   work happens. That is what makes the O(nActive * nHorizon) contract
            %   countable rather than merely reviewable (see the class header), and
            %   it is the point the anisotropic extension would have to change
            %  .
            %
            %   `rows` IS REQUIRED AND NEVER DEFAULTED. It is the whole active-set
            %   distinction made explicit at the one place it acts: `add` and
            %   `rebuild` pass every row because `dRaw` must stay true, `probe` and
            %   `probeGradient` pass the active rows because a resolved point
            %   cannot contribute. A default would hide which of those a call site
            %   meant, and getting it wrong in either direction is silent.
            %
            %   Protected rather than private so a test subclass can override it to
            %   count.
            %
            %   Inputs:
            %       Znorm  (nRows x nDim)  the samples to measure against
            %       rows   (nSup x 1) logical  which supporting points to measure
            D = rid.ext.pairwiseDistance(obj.Q(rows, :), Znorm);
        end

    end

    methods (Access = private)

        function d = nearestDistance(obj, Znorm)
            %NEARESTDISTANCE  Raw NN distance of EVERY supporting point, (nSup x 1).
            %
            %   All rows, always: this serves `add` and `rebuild`, which maintain
            %   the cache rather than exploit it.
            d = min(obj.distanceMatrix(Znorm, true(obj.nSup, 1)), [], 2);
        end

        function dRawCand = candidateDistance(obj, Zcand)
            %CANDIDATEDISTANCE  dRaw as it would stand with Zcand added, no mutation.
            %
            %   Measures only the active rows and leaves the rest at their cached
            %   value, which is exact: a resolved point's psi is zero either way.
            %   Shared by `probe` and nothing else -- `probeGradient` needs the
            %   argmin as well and does its own selection.
            dRawCand = obj.dRaw;
            act      = obj.isActive;
            if isempty(Zcand) || ~any(act)
                return
            end
            dRawCand(act) = min(obj.dRaw(act), ...
                                min(obj.distanceMatrix(Zcand, act), [], 2));
        end

        function refreshActive(obj)
            %REFRESHACTIVE  Recompute membership: active means NOT yet resolved.
            %
            %   Resolved is `dTilde_j <= 1`, i.e. `dRaw_j <= epsReq_j`. Exactly the
            %   complement, on the same non-strict boundary the criterion uses --
            %   a point at `dTilde = 1` has `psi = 0` and is not in the shell, so it
            %   is resolved.
            %
            %   Called by `add`, `setResolution` and `rebuild`; never by `probe`.
            obj.isActive = obj.dRaw > obj.epsReq;
        end

        function J = aggregate(obj, dRaw)
            %AGGREGATE  dRaw -> J.
            %
            %   Written out here rather than shared with rid.crit.evaluate, on
            %   purpose: a helper called by both would put the oracle and the
            %   subject under test on the same code path.
            dTilde = dRaw ./ obj.epsReq;
            psi    = min(max(dTilde - 1, 0), obj.psiMax);
            if isinf(obj.pNorm)
                J = max(psi);
            else
                J = mean(psi.^obj.pNorm)^(1 / obj.pNorm);
            end
        end

        function assertResolution(obj, epsReq)
            %ASSERTRESOLUTION  Shape and positivity of epsReq.
            if numel(epsReq) ~= obj.nSup
                error("rid:crit:invalidResolution", ...
                    "epsReq must have one entry per supporting point " + ...
                    "(got %d for nSup = %d).", numel(epsReq), obj.nSup);
            end
            if ~all(isfinite(epsReq)) || any(epsReq <= 0)
                error("rid:crit:invalidResolution", ...
                    "epsReq must be finite and strictly positive.");
            end
        end

        function assertColumns(obj, Znorm, name)
            %ASSERTCOLUMNS  Znorm lives in the same space as Q.
            if size(Znorm, 2) ~= obj.nDim
                error("rid:crit:dimensionMismatch", ...
                    "%s must have nDim = %d columns (got %d).", ...
                    name, obj.nDim, size(Znorm, 2));
            end
        end

    end

end
