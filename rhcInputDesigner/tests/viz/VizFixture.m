classdef VizFixture
%VIZFIXTURE  Runs at nDim = 2, 3 and 5, for the +viz tests.
%
%   THE FIGURES ARE ASSERTED AGAINST REAL RUNS, not against hand-built structs.
%   A +viz entry point takes a result or a record and reads only their declared
%   accessors, so a fabricated struct would test the drawing code against a
%   layout `rid.rhc.run` never produces -- which is the shape of a green
%   artefact that protects nothing.
%
%   TWO PROBLEMS, NOT ONE. `quiet` runs a plain budget design and produces no
%   events at all, so a progress plot of it must be ONE segment with NO marks.
%   `busy` is deliberately in the coverage reset's degenerate regime -- epsReq0 far above
%   what the horizon needs -- so it resets at nearly every step and refits its
%   surrogate on a schedule; it carries a REJECTED surrogate update as well as
%   applied ones, which is what separates "an event was logged" from "the J
%   baseline was broken".
%
%   The dimensions are reached through rid.rep.Narx alone: nDim = nInputs *
%   nLagU + nLagY, so (1,1) without the input difference is 2, (2,1) is 3 and
%   (3,2) is 5. The ROI is built from `requiredNames`, so nothing here has to
%   know the coordinate order.

    methods (Static)

        function r = quiet(nDim, extra)
            %QUIET  A plain offline budget run: a signal, and no events.
            arguments
                nDim  (1,1) double
                extra (1,1) struct = struct()
            end
            r = VizFixture.execute(nDim, 4, [], VizFixture.merge( ...
                struct("updateInterval", 0, ...
                       "epsReq0", VizFixture.threeStateBaseline(nDim)), extra));
        end

        function r = busy(extra)
            %BUSY  THE COVERAGE RESET's degenerate regime: resets, refits, and a rejected fit.
            arguments
                extra (1,1) struct = struct()
            end
            plant = FirstOrderProcess(0.8, 0.1, 1.4);
            r = VizFixture.execute(3, 8, plant, VizFixture.merge( ...
                struct("updateInterval", 2, "epsReq0", 1.2), extra));
        end

        function epsReq0 = threeStateBaseline(nDim)
            %THREESTATEBASELINE  A baseline at which ALL THREE states are occupied.
            %
            %   TUNED ON PURPOSE, and this is the fixture's whole job.
            %   At epsReq0 = 0.35 the two- and three-dimensional runs
            %   end with NO saturated point and the five-dimensional one with
            %   none resolved, so a two-colour resolved/unresolved plot agrees
            %   with a three-state one by construction and the state test cannot
            %   fail. These values were chosen by sweeping the baseline until
            %   every state is occupied at every dimension;
            %   `aFixtureInWhichAllThreeStatesAreOccupied` is what keeps that
            %   true rather than leaving it to a comment.
            %
            %   The value is a FIXTURE CHOICE and nothing else: it is not a
            %   recommended baseline, and the required epsReq0 leaves epsReq0 to the user.
            switch nDim
                case 2, epsReq0 = 0.18;
                case 3, epsReq0 = 0.18;
                case 5, epsReq0 = 0.40;
                otherwise
                    error("VizFixture:unsupportedDimension", ...
                        "nDim %d is not one of 2, 3, 5.", nDim);
            end
        end

        function r = saturating(extra)
            %SATURATING  A run the amplitude box binds at SOME samples, not all.
            %
            %   The fixture the saturation marks need, and both halves of that
            %   sentence are load-bearing. With no sample on a
            %   bound, marking and not marking agree; with EVERY sample on one,
            %   a mark-everything bug agrees too. At this baseline the
            %   two-dimensional design commits two of its four samples at
            %   u = -1 and two strictly inside, so only the correct rule passes.
            %
            %   The baseline is a FIXTURE CHOICE, not a recommendation: the required epsReq0
            %   leaves epsReq0 to the user.
            arguments
                extra (1,1) struct = struct()
            end
            r = VizFixture.execute(2, 4, [], VizFixture.merge( ...
                struct("updateInterval", 0, "epsReq0", 0.12), extra));
        end

        function r = saturatingWithPrimedBound()
            %SATURATINGWITHPRIMEDBOUND  ...and an INITIAL sample on the other bound.
            %
            %   `saturating` alone cannot show the priming exclusion to be
            %   wrong: its initial condition never touches a bound, so counting
            %   those samples and not counting them give the same number.
            %   Here the initial condition ends at u = +1
            %   while the design commits at u = -1, so the two are separable by
            %   sign as well as by count, and a rule that folds the initial
            %   condition into the design is off by exactly one.
            %
            %   `u(1,:)` is still the pre-run actuator position; it is
            %   the LAST priming row that sits on the bound, so the rate limit
            %   between the initial condition and the first designed sample is
            %   the only thing this changes about the problem.
            r = VizFixture.execute(2, 4, [], VizFixture.merge( ...
                struct("updateInterval", 0, "epsReq0", 0.12, ...
                       "uPrime", [0; 0.4; -0.3; 0.6; 1]), struct()));
        end

        function r = refitting(extra)
            %REFITTING  A run whose SURROGATE moves but whose coverage never resets.
            %
            %   The fixture the segment rule needs. In `busy` a coverage reset
            %   fires at nearly every step, so every segment is one step long
            %   whatever the surrogate rule says -- a rule keyed on the event
            %   KIND and a rule keyed on `detail.status` produce identical
            %   segments there, and the test cannot tell them apart.
            %   Here the only boundaries are applied surrogate
            %   updates, and the run also carries a REJECTED one, so the two
            %   rules disagree by construction.
            arguments
                extra (1,1) struct = struct()
            end
            plant = FirstOrderProcess(0.8, 0.1, 1.4);
            r = VizFixture.execute(3, 8, plant, VizFixture.merge( ...
                struct("updateInterval", 2, "epsReq0", 0.25), extra));
        end

        function rep = representation(nDim)
            %REPRESENTATION  The NARX space of the requested dimension.
            switch nDim
                case 2
                    rep = rid.rep.Narx(1, 1, "inputDifference", false, ...
                                       "pastSource", "simulated");
                case 3
                    rep = rid.rep.Narx(2, 1, "pastSource", "simulated");
                case 5
                    rep = rid.rep.Narx(3, 2, "pastSource", "simulated");
                otherwise
                    error("VizFixture:unsupportedDimension", ...
                        "nDim %d is not one of 2, 3, 5.", nDim);
            end
        end

        function roi = region(rep)
            %REGION  One declared range per required coordinate, by kind.
            roi = struct();
            for name = rep.requiredNames()
                if startsWith(name, "du")
                    roi.(name) = [-0.5 0.5];
                elseif startsWith(name, "y")
                    roi.(name) = [-2 2];
                else
                    roi.(name) = [-1 1];
                end
            end
        end

    end

    methods (Static, Access = private)

        function r = execute(nDim, nSteps, plant, cfg)
            %EXECUTE  Assemble a coherent problem and run it.
            epsReq0 = cfg.epsReq0;
            rep     = VizFixture.representation(nDim);
            roi = VizFixture.region(rep);

            [box, dimName] = rep.bounds(roi);
            nrm = rid.roi.Normalization(box, dimName);
            srg = rid.surrogate.LinearFirstOrder(0.5, 0.1, 1);

            nSup   = cfg.nSup;
            [Q, ~, qInfo] = rid.roi.supportingPoints(nrm, nSup, "seed", 1, ...
                                                      "progress", "off");
            fields = rid.fields.nominalFields(nSup);
            fillQ  = rid.roi.estimateFillDistance(Q, nrm);
            epsReq = rid.crit.demandedResolution(fields, "epsReq0", epsReq0, ...
                         "mode", "budget", "fillDistanceQ", fillQ);
            state  = rid.crit.CoverageState(Q, epsReq);

            gen = rid.gen.Direct(rid.gen.Constraints(rep, roi), 5);
            sol = rid.opt.MultiStart("seed", 1);
            if cfg.updateInterval > 0
                sch = rid.surrogate.ImmediateScheduler( ...
                          "updateInterval", cfg.updateInterval, "minSamples", 3);
            else
                sch = rid.surrogate.ImmediateScheduler();
            end
            setup = rid.rhc.Setup(srg, rep, nrm, roi, Q, fields, gen, sol, ...
                                  state, sch, rid.fields.ImmediateScheduler(), qInfo);

            ov = struct( ...
                "crit", struct("epsReq0", epsReq0), ...
                "rhc",  struct("nSteps", nSteps), ...
                "io",   struct("saveRun", false, "storeQ", cfg.storeQ));

            % Long enough that the first horizon is fully representable at every
            % dimension: the deepest space here discards three samples as warm-up.
            uPrime = cfg.uPrime;
            hist0  = struct("u", uPrime, "x0", srg.steadyState(uPrime(1)), ...
                            "tS", srg.tS);
            if ~isempty(plant)
                hist0.y = [0; 0.1; 0.05; 0.2; 0.1];
            end

            r = rid.rhc.run(plant, setup, hist0, ov);
        end

        function cfg = merge(base, extra)
            %MERGE  Defaults the fixtures share, overridden by a caller's keywords.
            cfg = base;
            if ~isfield(cfg, "storeQ")
                cfg.storeQ = true;
            end
            if ~isfield(cfg, "nSup")
                cfg.nSup = 8;
            end
            if ~isfield(cfg, "uPrime")
                cfg.uPrime = [0; 0.4; -0.3; 0.6; -0.5];
            end
            for name = string(fieldnames(extra)).'
                cfg.(name) = extra.(name);
            end
        end

    end

end
