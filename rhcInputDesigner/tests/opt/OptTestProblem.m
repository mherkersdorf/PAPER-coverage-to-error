classdef OptTestProblem < handle
%OPTTESTPROBLEM  A scoreFcn double that counts its own calls.
%
%   Contract item 1 of the method description -- `scoreFcn` is a black box
%   -- is not testable by inspection. It is testable by making the objective
%   observable in two ways an optimizer cannot fake:
%
%     - it COUNTS ITS OWN CALLS, so `diag.nEvaluations` must agree exactly. An
%       optimizer that evaluated anything outside the handle, or that cached and
%       skipped calls it claimed to make, disagrees.
%     - it carries a SECRET OFFSET added to every value. An optimizer that
%       reimplemented or short-circuited the criterion would return a best score
%       without it. No amount of knowledge about `J` can produce that number.
%
%   A handle class rather than a closure, so the count survives being passed
%   around and can be read after the search.
%
%   The three arities of the method description are separate methods on
%   purpose. `withGradient` can be asked for one, two or three outputs, so it
%   exercises the analytic path; `valueOnly` can be asked for exactly one, so an
%   optimizer must fall back to differencing it -- which is what makes the
%   gradient-free path a real test rather than a configuration flag.
%
%   `nShellFcn` supplies `info.nShell`, the count of supporting points in the shell
%   at that theta. It is a function of theta and not a constant, because that is
%   what it is in the real objective: the dead-on-arrival diagnostic asks whether
%   THIS start sits in a flat region, and a constant would make every start alike.

    properties
        count  (1,1) double = 0
        offset (1,1) double = 0        % the secret; see the header
        objective                      % theta -> f, the true objective
        gradient                       % theta -> g, its gradient
        nShellFcn                      % theta -> nShell, supporting points in shell
    end

    methods

        function obj = OptTestProblem(objective, gradient, offset, nShellFcn)
            arguments
                objective (1,1) function_handle
                gradient  (1,1) function_handle
                offset    (1,1) double = 0
                % The default says "some point is in the shell", so a problem that
                % is not about dead-on-arrival never reports a dead start by
                % accident.
                nShellFcn (1,1) function_handle = @(theta) 1
            end
            obj.objective = objective;
            obj.gradient  = gradient;
            obj.offset    = offset;
            obj.nShellFcn = nShellFcn;
        end

        function [f, g, info] = withGradient(obj, theta)
            %WITHGRADIENT  The [f], [f, g] and [f, g, info] forms.
            obj.count = obj.count + 1;
            f = obj.objective(theta) + obj.offset;
            if nargout > 1
                g = obj.gradient(theta);
            end
            if nargout > 2
                % J, the penalty and nShell reported SEPARATELY.
                % Here the offset plays the part of the penalty: it is the
                % difference between the criterion and the objective, which is
                % exactly what a penalty is.
                info = struct("J", f - obj.offset, "penalty", obj.offset, ...
                              "nShell", obj.nShellFcn(theta));
            end
        end

        function f = valueOnly(obj, theta)
            %VALUEONLY  The gradient-free form: one output and no more.
            obj.count = obj.count + 1;
            f = obj.objective(theta) + obj.offset;
        end

    end

    methods (Static)

        function gen = generator(nHorizon, amplitude, rate)
            %GENERATOR  A single-input rid.gen.Direct over a declared box.
            %
            %   The REAL generator, not a double: `solve` is specified against
            %   `feasibleSet`, `toTrajectory`, `nParams` and `initialGuesses`, and
            %   a stub of those would test the stub.
            arguments
                nHorizon  (1,1) double
                amplitude (1,2) double
                rate      (1,2) double
            end
            rep = rid.rep.Narx(2, 1);
            roi = struct("u1_k", amplitude, "du1_k", rate, "y_km1", [-2 3]);
            gen = rid.gen.Direct(rid.gen.Constraints(rep, roi), nHorizon);
        end

        function fcn = fixedStarts(Theta0)
            %FIXEDSTARTS  A starts provider returning exactly these rows.
            %
            %   The dead-on-arrival test needs a KNOWN number of starts in the flat
            %   region, which a random provider cannot give. It honours the
            %   provider contract -- nStarts in, [Theta0, info] out -- and refuses
            %   a count it was not built for rather than silently returning the
            %   wrong number of rows.
            fcn = @provide;
            function [T, info] = provide(nStarts)
                if nStarts ~= size(Theta0, 1)
                    error("rid:opt:invalidStarts", ...
                        "This provider holds %d start(s); %d were asked for.", ...
                        size(Theta0, 1), nStarts);
                end
                T    = Theta0;
                info = struct("seed", NaN, "nStarts", nStarts, ...
                              "warmStart", "none", "source", "fixed");
            end
        end

    end

end
