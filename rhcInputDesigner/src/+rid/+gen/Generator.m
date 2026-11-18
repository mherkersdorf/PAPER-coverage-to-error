classdef Generator
%GENERATOR  The contract every input parametrization must satisfy.
%
%   A parametrization decides what a candidate IS -- which numbers the optimizer
%   moves, what trajectory they denote, and how much of that trajectory is
%   committed per iteration. The loop and the optimizer talk to this interface and
%   never to a particular parametrization:
%
%       gen = rid.gen.Direct(constraints, nHorizon)   % every sample a variable
%       gen = rid.gen.BSpline(...)                    % a reduced one
%
%   THERE IS NO CANDIDATE SET, and that is what this interface exists to express.
%   For the default parametrization the decision space is `nHorizon * nInputs`
%   dimensional, so a generator hands out a DESCRIPTION of the feasible
%   set and a handful of feasible starting points -- never a list.
%
%   WHY IT EXISTS AT ALL, since `rid.gen.Direct` is the only implementation. Two
%   reasons, and neither is symmetry for its own sake:
%
%     - `rid.rhc.Setup` validates its nine arguments BY CONCRETE CLASS.
%       Without a base class the generator was the one argument that could not be
%       typed, and it is the one most easily swapped for another object, because
%       it sits between `Q` and `state` in a nine-argument positional signature.
%       Duck-typing would make the contract a list inside an error message rather
%       than a declaration the language enforces.
%     - the reduced parametrizations need it regardless, and building it
%       around one implementation is cheaper than retrofitting it around two.
%
%   Every other component in the project is already this pattern --
%   `rid.surrogate.Surrogate`, `rid.rep.Representation`,
%   `rid.surrogate.Scheduler`, `rid.fields.Scheduler`, `rid.opt.Optimizer`. The
%   generator being the exception was an accident of build order, not a design.
%
%   THE APPLIED UNIT IS DECLARED HERE, NOT ASSUMED BY THE LOOP. `applied`
%   returns one sample for `rid.gen.Direct` and the first segment for a spline;
%   `+rhc` commits exactly what it returns and never "the first sample". That one
%   method is what lets a new parametrization enter without touching the loop.
%
%   PHYSICAL UNITS throughout. A generator works in the physical units of the
%   actuator and of the declared admissible region; normalization happens at the
%   representation boundary, in the loop. `penalty` is the one place a generator
%   sees regressor coordinates, and they are physical there too.
%
%   GENERATORS ARE VALUE OBJECTS. They hold a declaration and no run state, so the
%   two intended homes of mutable state stay `rid.crit.CoverageState`
%   and the update schedulers.
%
%   NO CONSTRUCTOR IS DECLARED HERE. What a parametrization needs in order to
%   exist differs between parametrizations -- `rid.gen.Direct` takes a constraint
%   declaration and a horizon, a spline generator would take a knot vector as well
%   -- and the base class has nothing of its own to store. Contrast
%   `rid.rep.Representation`, which owns `pastSource` and therefore does have a
%   constructor every subclass forwards to.

    properties (Abstract, SetAccess = immutable)
        % Declared without size/class validators: MATLAB does not permit them on a
        % property that implements an inherited abstract one, and repeating them in
        % every subclass is not an option either. Shape:
        %
        %   nHorizon  (1x1) double  prediction horizon in samples
        %
        % USER-DEFINED AND FIXED FOR THE RUN, with
        % no default anywhere. The generator is its SINGLE HOLDER: the
        % surrogate's simulation length, the representation's future coordinates
        % and the penalty's `gamma` all read this one value, and `+rhc` reads it
        % here rather than keeping a copy.
        %
        % It is a property rather than one of the abstract methods below because
        % the loop reads it once per step to size the candidate window and
        % a method call would suggest it could change.
        nHorizon

        %   constraints  (1x1) rid.gen.Constraints  the hard amplitude and rate
        %                                           box, PHYSICAL
        %
        % ON THE CONTRACT BECAUSE THE LIMITS BELONG TO THE ACTUATOR, not
        % to the parametrization. Every parametrization that designs an input
        % signal is subject to the same physical box -- the reduced ones of
        % the method description included. What IS
        % parametrization-specific is how the box is enforced in a given
        % parameter space, and that is exactly what `feasibleSet` returns; the
        % two are different questions and only the second varies.
        %
        % WHY IT IS NOT LEFT TO EACH SUBCLASS. A generator without the property
        % does not raise: a consumer asking for the box gets nothing and
        % continues. `rid.viz.signal` then draws no bounds, and a run pinned
        % against its actuator limits becomes indistinguishable from one that
        % never approached them -- a silent absence that changes what a reader
        % concludes, which is the failure mode the generator contract introduced this class to
        % remove.
        constraints
    end

    methods (Abstract)

        % p = nInputs(obj)
        %   p  (1x1) double  input channels. A METHOD, not a property: the channel
        %                    count has one source -- the representation, through
        %                    the constraint declaration -- and a stored copy could
        %                    disagree with it.
        p = nInputs(obj)

        % n = nParams(obj)
        %   n  (1x1) double  length of the decision vector. Takes no argument
        %                    beyond obj: every quantity it needs lives on the
        %                    object already.
        n = nParams(obj)

        % uSeq = toTrajectory(obj, theta)
        %   theta  (nParams x 1)        the decision vector
        %   uSeq   (nHorizon x nInputs) PHYSICAL input trajectory
        uSeq = toTrajectory(obj, theta)

        % uApp = applied(obj, uSeq)
        %   uSeq  (nHorizon x nInputs)  a full candidate trajectory
        %   uApp  (nApplied x nInputs)  THE UNIT COMMITTED THIS ITERATION.
        %                               One sample for rid.gen.Direct; the first
        %                               segment for a spline. The loop commits
        %                               exactly this and never assumes its height.
        uApp = applied(obj, uSeq)

        % theta = shift(obj, thetaPrev)
        %   The warm start: the previous solution advanced by one applied
        %   unit. FEASIBLE BY CONSTRUCTION -- there is no repair path anywhere in
        %   this package, and a start that violates the feasible set means the
        %   construction is wrong, not the start.
        %
        %   thetaPrev  (nParams x 1)
        %   theta      (nParams x 1)
        theta = shift(obj, thetaPrev)

        % [lb, ub, A, b] = feasibleSet(obj, hist)
        %   The HARD constraints, exactly: bounds and linear inequalities, so an
        %   infeasible candidate is unrepresentable rather than penalized
        %   (the project's coding standard invariant 6). STATE-DEPENDENT through u(k-1), which
        %   is why it takes the history.
        %
        %   hist  (1x1) struct           the signal history
        %   lb    (nParams x 1) double
        %   ub    (nParams x 1) double
        %   A     (nRows x nParams)      may be sparse; empty when there are none
        %   b     (nRows x 1) double
        [lb, ub, A, b] = feasibleSet(obj, hist)

        % [Theta0, info] = initialGuesses(obj, hist, thetaPrev, opts)
        %   The multi-start set, every row feasible by construction. `nStarts` is a
        %   REQUIRED name-value argument with no default here: how many local
        %   searches to run is the optimizer's decision and is declared at
        %   `opts.opt.nStarts` -- +opt decides the count, +gen answers.
        %
        %   hist       (1x1) struct
        %   thetaPrev  (nParams x 1)     previous solution; EMPTY on the first step
        %   Theta0     (nStarts x nParams)
        %   info       (1x1) struct      what the construction did; at minimum
        %                                .seed .nStarts .warmStart
        [Theta0, info] = initialGuesses(obj, hist, thetaPrev, opts)

        % [pen, gradPen, info] = penalty(obj, Zphys, opts)
        %   The soft-constraint term of the OBJECTIVE, never a term in `J`.
        %   `pen` is the WEIGHTED total, so the caller adds it directly, and
        %   `gradPen` is with respect to `Zphys` rather than to `theta` so that it
        %   shares the simulation Jacobian with the criterion's gradient path.
        %
        %   Zphys    (nHorizon x nDim)  PHYSICAL regressor coordinates
        %   pen      (1x1) double
        %   gradPen  (nHorizon x nDim)  d(pen)/d(Zphys)
        %   info     (1x1) struct
        [pen, gradPen, info] = penalty(obj, Zphys, opts)

        % s = describe(obj)
        %   s  (1x1) string  human-readable summary, for logs and run records
        s = describe(obj)

    end

end
