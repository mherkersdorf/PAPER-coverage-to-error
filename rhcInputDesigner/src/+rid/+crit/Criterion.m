classdef (Abstract) Criterion < handle
%CRITERION  What the receding-horizon loop asks of a criterion.
%
%   The loop is handed a criterion and calls it `state`. Every
%   implementation scores POINTS -- normalized regressor-space rows, samples in
%   rows -- and never sees a model: `+crit` must not depend on `+surrogate` or
%   `+rep`. The criterion component interface's `ctx` argument is therefore not
%   taken; if a model-based criterion ever arrives, the retrofit touches
%   `stepObjective` and `rid.rhc.run` only.
%
%   SEVEN ABSTRACT METHODS, which every criterion must define:
%
%       add(Znew)                         commit realized points
%       J = probe(Zcand)                  score a candidate; NEVER mutates
%       [gradCand, info] = probeGradient(Zcand)
%                                         dJ/dZcand, (nCand x nDim); info.J and
%                                         info.nShell
%       J = value()                       the criterion at the committed points
%       rebuild(Zall)                     recompute from these points; empty
%                                         means nothing committed
%       fresh = freshCopy()               same configuration, nothing committed
%       s = describe()                    one line for the run record
%
%   Scores are MINIMIZED; an adapter negates a criterion that is
%   naturally maximized.
%
%   `info.nShell` IS THE NUMBER OF TERMS THE CRITERION GRADIENT IS A SUM OVER, so
%   `nShell == 0` holds exactly when the criterion gradient is exactly zero. It is
%   always a count, never NaN: the optimizer reads dead-on-arrival from it,
%   and `NaN == 0` would report a start as
%   alive without anyone having asked. For the coverage criterion it is the
%   number of shell points; every implementation states its own.
%
%   SEVEN DEFAULTS, each the honest answer of a criterion that LACKS the concept.
%   They are how the loop reaches coverage-specific features without asking
%   which class it holds -- no `isa` anywhere in the loop:
%
%       isSatisfied()        false   certification is a property of a criterion,
%                                    not of criteria in general. False
%                                    means: never "certified", and the budget-mode
%                                    reset of the coverage reset never fires
%       isMonotone()         false   `value()` is non-increasing under `add`
%                                    between rebuilds and resolution changes;
%                                    gates invariant 3
%       resolution()         empty   the demanded resolution `epsReq`; empty means
%                                    the criterion has none, and the record then
%                                    stores NaN
%       setResolution(e)     raises  `rid:crit:noResolution`: a field update for
%                                    a criterion without a resolution is an error,
%                                    never silently dropped
%       nActive()            NaN     supporting points still in play
%       isConsistentWith(Q)  true    whether the criterion was built against
%                                    these supporting points (`rid.rhc.Setup`)
%       withoutQ()           obj     a copy with the supporting points removed
%                                   ; nothing to remove here

    methods (Abstract)
        add(obj, Znew)
        J = probe(obj, Zcand)
        [gradCand, info] = probeGradient(obj, Zcand)
        J = value(obj)
        rebuild(obj, Zall)
        fresh = freshCopy(obj)
        s = describe(obj)
    end

    methods

        function tf = isSatisfied(obj) %#ok<MANU>
            %ISSATISFIED  Has this criterion certified its own completion? No.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = false;
        end

        function tf = isMonotone(obj) %#ok<MANU>
            %ISMONOTONE  Is value() non-increasing under add? Not unless declared.
            %
            %   Outputs:
            %       tf  (1x1) logical
            tf = false;
        end

        function epsReq = resolution(obj) %#ok<MANU>
            %RESOLUTION  The demanded resolution; empty for a criterion with none.
            %
            %   Outputs:
            %       epsReq  (0 x 1) double
            epsReq = zeros(0, 1);
        end

        function setResolution(obj, epsReq) %#ok<INUSD>
            %SETRESOLUTION  Refused: this criterion has no demanded resolution.
            error("rid:crit:noResolution", ...
                "%s has no demanded resolution, so a refinement-field update " + ...
                "cannot act on it. A field scheduler that fires for such a " + ...
                "criterion is a configuration error, not something to drop " + ...
                "silently.", class(obj));
        end

        function n = nActive(obj) %#ok<MANU>
            %NACTIVE  Supporting points still in play; NaN, the concept is absent.
            %
            %   Outputs:
            %       n  (1x1) double
            n = NaN;
        end

        function tf = isConsistentWith(obj, Q) %#ok<INUSD>
            %ISCONSISTENTWITH  Built against these supporting points? Nothing to check.
            %
            %   Inputs:
            %       Q   (nSup x nDim) double
            %   Outputs:
            %       tf  (1x1) logical
            tf = true;
        end

        function copy = withoutQ(obj)
            %WITHOUTQ  The criterion with its supporting points removed: it has none.
            %
            %   Outputs:
            %       copy  (1x1) rid.crit.Criterion  this object
            copy = obj;
        end

    end

end
