classdef Setup
%SETUP  Everything that defines the PROBLEM, as one object.
%
%       setup  = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, gen, sol, state, ...
%                              sch, fieldSch, qInfo)
%       result = rid.rhc.run(process, setup, hist0, opts)
%
%   THE GROUPING IS NOT TIDINESS. `setup` is exactly what a run record must store
%   as "what problem was this" and exactly what
%   `replay` consumes. When the grouping the code wants and the grouping the record
%   wants coincide, that is the grouping. The alternative was a ten-argument
%   positional signature that the refinement-field scheduler would have made eleven.
%
%   WHAT STAYS OUTSIDE, AND WHY. `process` is the one thing a replay must NOT have
%   -- a replay re-derives from the stored history and must never touch hardware --
%   so keeping it out of `Setup` makes that structural rather than remembered.
%   `hist0` stays outside because it is the experiment's INITIAL CONDITION: data,
%   in the plain history struct, and not part of the problem definition.
%   `opts` stays outside because it is values rather than objects, and the record
%   already stores it separately as an effective struct plus the overrides.
%
%   WHAT CAME IN, AND WHY. `sol` is here because **a replay must
%   have the same solver or it does not reproduce** -- the exact mirror of the
%   argument that keeps `process` out. Constructed inside `run` from `opts` it
%   made `rid.opt.GridSearch` unselectable, and the method description's
%   "checkable against exhaustive search" true of `+opt` in isolation and never of
%   a run. `fields` is here because three `(nSup x 1)` vectors are DATA living
%   beside `Q`, not settings; routing them through `opts` contradicted this
%   entry's own "values, not objects".
%
%   THE CRITERION IS A HANDLE, AND THAT HAS A RULE. `state` is a
%   `rid.crit.Criterion`, which the loop MUTATES, so storing `setup` at
%   the end of a run does not capture the configuration it started from: the
%   record stores `epsReqInitial` beside it, NaN for
%   a criterion with no resolution, whose configuration nothing in the loop
%   moves. The rule generalises: any handle member arriving here needs the same
%   treatment.
%
%   NOT NAME-VALUE. Name-value is order-proof and self-documenting, but wordy at
%   every call site, and a mistyped name silently drops an argument instead of
%   failing. The `arguments` block below validates all nine BY CONCRETE CLASS,
%   which makes almost every misordering an immediate typed error -- and is why
%   `rid.gen.Generator` had to exist first: without it the ninth argument
%   was the one that could not be typed, and it sits between `Q` and `state`.
%
%   A TWELFTH ARGUMENT, `qInfo`, AND IT HAS NO DEFAULT (the Q provenance, amending the single-setup rule).
%   Where `Q` came from is part of "what problem was this", and an absent
%   provenance is precisely the state this argument removes: a `Setup` built
%   without it does not construct. `"supplied"` is a legitimate value the caller
%   passes EXPLICITLY when `Q` did not come from `rid.roi.supportingPoints` --
%   never a fallback the code invents when the caller omits the argument.
%
%   THE COHERENCE CHECKS ARE THE OTHER HALF. Typing catches a misordering; it does
%   not catch a `Q` built for a different normalization, or a coverage state
%   constructed against different supporting points. Both are silent: the run
%   proceeds and measures coverage of a space nobody declared. They are checked
%   here, once, at the only moment the objects are together and before anything has
%   been computed. `rep.assertCompatible(srg)` is the same argument across the
%   representation/surrogate pair: a space whose coordinates ARE the surrogate's
%   states is incoherent with a surrogate of another order, and only the
%   representation knows whether it has such a precondition.
%
%   WHAT IS NOT CHECKED HERE. The `process = []` with `pastSource = "measured"`
%   refusal belongs to `rid.rhc.run`, because
%   `process` is deliberately not in this object. Putting the check here would
%   require the argument the design keeps out.
%
%   Construction:
%       setup = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, gen, sol, state, ...
%                             sch, fieldSch, qInfo)
%
%   Inputs:
%       srg       (1x1) rid.surrogate.Surrogate   the geometric guide
%       rep       (1x1) rid.rep.Representation    the coordinates coverage is in
%       nrm       (1x1) rid.roi.Normalization     physical <-> normalized
%       roi       (1x1) struct                    one (1 x 2) range per dimension
%       Q         (nSup x nDim) double            supporting points, NORMALIZED
%       fields    (1x1) struct                    .relevance .complexity
%                                                 .noiseLevel, each (nSup x 1)
%       gen       (1x1) rid.gen.Generator         the parametrization
%       sol       (1x1) rid.opt.Optimizer         the search
%       state     (1x1) rid.crit.Criterion        the criterion, named `state`
%                                                
%       sch       (1x1) rid.surrogate.Scheduler   surrogate update policy
%       fieldSch  (1x1) rid.fields.Scheduler      field update policy
%       qInfo     (1x1) struct                    Q's provenance, fields
%                                                 .origin ("generated" | "supplied"
%                                                 | "stored"), .method, .seed,
%                                                 .hash -- "" / NaN where the
%                                                 concept does not apply. The
%                                                 struct rid.roi.supportingPoints
%                                                 returns as its third output,
%                                                 tagged with .origin; built by
%                                                 hand for a Q that never went
%                                                 through it.

    properties (SetAccess = immutable)
        % Declared WITHOUT class validators, and the typing lives in the
        % `arguments` block below -- which is what the single-setup rule specifies and, as it
        % happens, the only thing MATLAB permits: a property validated against an
        % abstract class needs a default value it cannot have, since the class
        % cannot be instantiated to produce one. Nothing is lost. The constructor
        % is the only writer, the properties are immutable, and every argument is
        % typed on the way in.
        %
        %   srg       (1x1) rid.surrogate.Surrogate
        %   rep       (1x1) rid.rep.Representation
        %   nrm       (1x1) rid.roi.Normalization
        %   roi       (1x1) struct
        %   Q         (nSup x nDim) double, NORMALIZED
        %   fields    (1x1) struct, three (nSup x 1) members
        %   gen       (1x1) rid.gen.Generator
        %   sol       (1x1) rid.opt.Optimizer
        %   state     (1x1) rid.crit.Criterion       -- a HANDLE; the loop advances it
        %   sch       (1x1) rid.surrogate.Scheduler  -- a HANDLE
        %   fieldSch  (1x1) rid.fields.Scheduler     -- a HANDLE
        %   qInfo     (1x1) struct  Q's provenance: .origin, .method,
        %                           .seed, .hash
        srg
        rep
        nrm
        roi
        fields
        gen
        sol
        sch
        fieldSch
        qInfo
    end

    properties (SetAccess = private)
        % THE TWO MEMBERS `withoutQ` HAS TO REACH, and private rather than
        % immutable for that reason alone.
        % `opts.io.storeQ = false` must remove `Q` from a stored record
        % completely, and the single-setup rule's "the setup is stored in full" kept two further
        % copies of it -- here and inside `state` -- so a record could announce
        % `storedQ = false` and carry the points regardless. Nothing else writes
        % either after construction; the constructor remains the only other
        % writer and the coherence checks still run there.
        %
        %   Q      (nSup x nDim) double            supporting points, NORMALIZED
        %   state  (1x1) rid.crit.Criterion        a HANDLE; the loop advances it
        Q
        state
    end

    methods

        function copy = withoutQ(obj)
            %WITHOUTQ  A copy of this setup with the supporting points removed.
            %
            %   THE ONE EXCEPTION TO "STORED IN FULL",
            %   and it is written down rather than left as a back door. The flag
            %   exists for SIZE: while `Q` sat inside the stored `Setup` the file
            %   never actually got smaller, so the option did nothing while
            %   announcing that it had.
            %
            %   A VALUE CLASS, so the copy is free -- but `state` is a handle,
            %   and rebinding it here points the COPY at a new object while
            %   leaving the caller's own state alone. Clearing in place would
            %   empty the run's live `result` as a side effect of recording it.
            %
            %   The consequence is the documented one: such a record is not
            %   exactly replayable, and `rid.io.replay` says so and stops.
            %
            %   Outputs:
            %       copy  (1x1) rid.rhc.Setup  with `Q` (0 x nDim) and a state
            %                                  that has none either
            copy       = obj;
            copy.Q     = zeros(0, size(obj.Q, 2));
            copy.state = obj.state.withoutQ();
        end

    end

    methods

        function obj = Setup(srg, rep, nrm, roi, Q, fields, gen, sol, state, ...
                             sch, fieldSch, qInfo)
            arguments
                srg      (1,1) rid.surrogate.Surrogate
                rep      (1,1) rid.rep.Representation
                nrm      (1,1) rid.roi.Normalization
                roi      (1,1) struct
                Q        (:,:) double {mustBeReal, mustBeFinite}
                fields   (1,1) struct
                gen      (1,1) rid.gen.Generator
                sol      (1,1) rid.opt.Optimizer
                state    (1,1) rid.crit.Criterion
                sch      (1,1) rid.surrogate.Scheduler
                fieldSch (1,1) rid.fields.Scheduler
                qInfo    (1,1) struct
            end

            assertValidQInfo(qInfo);

            % The ROI declaration is validated by the representation that owns the
            % coordinate names, which also gives the missing-range and
            % unknown-range diagnostics rather than a generic one here.
            [box, dimName] = rep.bounds(roi);

            if rep.nDim ~= nrm.nDim
                error("rid:rhc:inconsistentSetup", ...
                    "The representation has %d dimension(s) and the " + ...
                    "normalization %d. They describe the same space.", ...
                    rep.nDim, nrm.nDim);
            end
            % The normalization IS the ROI expressed as a metric, so a
            % normalization built from a different declaration silently measures
            % coverage in the wrong units -- distances scale, `epsReq` does not, and
            % the run reports a resolution it never achieved.
            if ~isequal(dimName(:), nrm.dimName(:))
                error("rid:rhc:inconsistentSetup", ...
                    "The normalization is over %s and the representation " + ...
                    "declares %s.", strjoin("""" + nrm.dimName + """", ", "), ...
                    strjoin("""" + dimName + """", ", "));
            end
            if ~isequal(box, nrm.box)
                bad = find(any(box ~= nrm.box, 2), 1);
                error("rid:rhc:inconsistentSetup", ...
                    "The ROI declares [%g, %g] for ""%s"" and the normalization " + ...
                    "was built over [%g, %g]. The normalization is the ROI " + ...
                    "expressed as a metric; two declarations of one box " + ...
                    "measure coverage in different units without raising.", ...
                    box(bad, 1), box(bad, 2), dimName(bad), ...
                    nrm.box(bad, 1), nrm.box(bad, 2));
            end

            % THE SPACE CAN BE COHERENT AND STILL NOT BE THIS SURROGATE'S. The
            % three checks above compare declarations with declarations; this one
            % asks the representation whether the model it will be built from can
            % supply its coordinates at all. It is a no-op for
            % every representation whose coordinates are arithmetic on recorded
            % signals, and rid.rep.StateSpace refuses a state dimension that is
            % not its own -- HERE, at assembly, rather than inside the first
            % `build` of a run that has already started a process. The
            % representation owns the diagnostic, so nothing here names a
            % representation class.
            rep.assertCompatible(srg);

            if size(Q, 2) ~= nrm.nDim
                error("rid:rhc:inconsistentSetup", ...
                    "Q has %d column(s); the space has %d dimension(s). " + ...
                    "Supporting points are NORMALIZED coordinates, samples in rows.", ...
                    size(Q, 2), nrm.nDim);
            end
            % A coverage state caches distances to the Q it was CONSTRUCTED with.
            % A different Q here would leave the loop reporting coverage of one set
            % while `resolutionReport` and the record described another. The
            % criterion answers; one without supporting points has nothing to
            % disagree with.
            if ~state.isConsistentWith(Q)
                error("rid:rhc:inconsistentSetup", ...
                    "The criterion was constructed against a different Q " + ...
                    "(%s here). It caches distances to its own supporting " + ...
                    "points, so the two cannot differ.", mat2str(size(Q)));
            end

            % Completeness is checked by the package that CONSTRUCTS fields, so
            % there is one definition of a complete set. Values --
            % finiteness and strict positivity -- belong to
            % rid.crit.demandedResolution, which owns them.
            rid.fields.Scheduler.assertCompleteFields(fields);
            for name = ["relevance", "complexity", "noiseLevel"]
                value = fields.(name);
                if ~isequal(size(value), [size(Q, 1), 1])
                    error("rid:rhc:inconsistentSetup", ...
                        "fields.%s is %s; it carries one value per supporting " + ...
                        "point, so it must be (%d x 1). The fields live beside " + ...
                        "Q and are indexed by it.", ...
                        name, mat2str(size(value)), size(Q, 1));
                end
            end

            obj.srg      = srg;
            obj.rep      = rep;
            obj.nrm      = nrm;
            obj.roi      = roi;
            obj.Q        = Q;
            obj.fields   = fields;
            obj.gen      = gen;
            obj.sol      = sol;
            obj.state    = state;
            obj.sch      = sch;
            obj.fieldSch = fieldSch;
            obj.qInfo    = qInfo;
        end

        function s = describe(obj)
            %DESCRIBE  One block for the run record and the console of a caller.
            %
            %   Outputs:
            %       s  (1x1) string
            s = "rid.rhc.Setup" + newline + ...
                "  surrogate:      " + obj.srg.describe() + newline + ...
                "  representation: " + obj.rep.describe() + newline + ...
                "  normalization:  " + obj.nrm.describe() + newline + ...
                "  generator:      " + obj.gen.describe() + newline + ...
                "  optimizer:      " + class(obj.sol) + newline + ...
                "  supporting:     " + sprintf("%d point(s) in %d dimension(s)", ...
                    size(obj.Q, 1), size(obj.Q, 2)) + newline + ...
                "  schedulers:     " + obj.sch.describe() + " / " + ...
                    obj.fieldSch.describe();
        end

    end

end

% ---------------------------------------------------------------------------

function assertValidQInfo(qInfo)
%ASSERTVALIDQINFO  Q's provenance is complete and, for a supplied Q, invents
%nothing.
%
%   `.origin` is CLOSED: "generated" (rid.roi.supportingPoints ran), "supplied"
%   (a hand-built Q the caller declares), "stored" (the stored-Q reuse, a previously
%   generated set loaded from disk). `.method` and `.seed` describe the
%   generator and are `"" ` / `NaN` where the concept does not apply; `.hash`
%   is the content hash of a stored set and `""` otherwise. A record
%   that omits is better than one that misreports, so "supplied" is
%   refused unless every other field is the explicit absence, never left for
%   the caller to forget.
required = ["origin", "method", "seed", "hash"];
missing  = required(~isfield(qInfo, required));
if ~isempty(missing)
    error("rid:rhc:invalidQInfo", ...
        "qInfo is missing %s. It must declare origin, method, seed and hash " + ...
        ": the run record reads Q's provenance from here and " + ...
        "nowhere else.", strjoin(missing, ", "));
end

origin = string(qInfo.origin);
if ~ismember(origin, ["generated", "supplied", "stored"])
    error("rid:rhc:invalidQInfo", ...
        "qInfo.origin must be ""generated"", ""supplied"" or ""stored"" " + ...
        "(got ""%s"").", origin);
end

if origin == "supplied" && ...
        (string(qInfo.method) ~= "" || ~isnan(qInfo.seed) || string(qInfo.hash) ~= "")
    error("rid:rhc:invalidQInfo", ...
        "qInfo.origin is ""supplied"", so .method and .hash must be """" and " + ...
        ".seed must be NaN: nothing about a supplied Q's provenance may be " + ...
        "invented.");
end

end
