classdef GeneratorContractTest < matlab.unittest.TestCase
%GENERATORCONTRACTTEST  The parametrization contract of rid.gen.Generator.
%
%   The base class exists so that `rid.rhc.Setup` can type its ninth argument by
%   concrete class rather than by a list inside an error message. What
%   makes that worth testing is not that `rid.gen.Direct` inherits -- a one-word
%   change -- but that the DECLARATION is complete: a class satisfying the whole
%   documented contract must be instantiable, and one missing any single member
%   must not be.
%
%   Each test names the input that would make it fail. A test that only asserted
%   `isa(rid.gen.Direct(...), "rid.gen.Generator")` would pass against a base
%   class declaring nothing at all, which is the failure mode worth catching:
%   an empty abstract block types the argument and enforces no contract.

    methods (Test)

        function directIsAGeneratorAndTheContractIsComplete(tc)
            % Fails if: Direct stops inheriting, or the abstract block names a
            % member Direct does not implement -- MATLAB would then refuse to
            % instantiate it, which is the enforcement the base class is for.
            gen = GeneratorContractTest.generator(5);
            tc.verifyClass(gen, "rid.gen.Direct");
            tc.verifyTrue(isa(gen, "rid.gen.Generator"));
        end

        function everyDocumentedMemberIsDeclaredAbstract(tc)
            % The generator contract enumerates the contract, and `rid.rhc.Setup` and the loop rely
            % on it. If a member is dropped from the abstract block the class still
            % instantiates, `Setup` still types the argument, and a second
            % parametrization can then omit the member entirely -- discovered at
            % run time, inside the loop.
            %
            % `nHorizon` is a PROPERTY rather than a method: the loop reads it once
            % per step to size the candidate window and it is fixed for the
            % run.
            %
            % Fails if: a member is removed from the abstract block, or nHorizon is
            % demoted to a concrete property of Direct alone.
            mc       = meta.class.fromName("rid.gen.Generator");
            abstractMethods = string({mc.MethodList([mc.MethodList.Abstract]).Name});
            expected = ["nInputs", "nParams", "toTrajectory", "applied", ...
                        "shift", "feasibleSet", "initialGuesses", "penalty", ...
                        "describe"];
            tc.verifyEmpty(setdiff(expected, abstractMethods), ...
                "a documented member is missing from the abstract block");

            % `constraints` joined the contract at the constraints contract: the amplitude and rate
            % box belongs to the ACTUATOR, and every parametrization that designs
            % an input signal is subject to the same physical limits. How they are
            % enforced in a given parameter space is what `feasibleSet` returns,
            % and that is the parametrization-specific half.
            abstractProps = string({mc.PropertyList([mc.PropertyList.Abstract]).Name});
            tc.verifyEqual(sort(abstractProps), ["constraints", "nHorizon"]);
        end

        function aGeneratorWithoutConstraintsCannotBeInstantiated(tc)
            % THE SHARP ONE FOR THE CONSTRAINTS CONTRACT, and it needs its own double. A generator
            % lacking the box used not to raise at all: the plot simply drew no
            % bounds, so a run pinned against its actuator limits looked exactly
            % like one that never approached them. A silent absence that changes
            % what a reader concludes is worse than an error.
            %
            % Fails if: `constraints` is removed from the abstract block, or
            % demoted back to a concrete property of rid.gen.Direct alone.
            tc.verifyError(@() ConstraintlessGeneratorDouble(), ...
                "MATLAB:class:abstract");
        end

        function theRefusalIsAboutConstraintsAndNothingElse(tc)
            % Anchors the assertion where the statement is sharp.
            % `MATLAB:class:abstract` fires for the nine methods too, and the
            % error MESSAGE names no member at all, so neither shows that the
            % PROPERTY block did the refusing. `abstractDetails` is what does:
            % it returns exactly the unimplemented members.
            missing = matlab.metadata.abstractDetails( ...
                "ConstraintlessGeneratorDouble");
            tc.verifyEqual(string({missing.Name}), "constraints", ...
                "the double implements everything else, so `constraints` must " + ...
                "be the only thing standing between it and instantiation");
            tc.verifyClass(missing, "matlab.metadata.Property", ...
                "and it must be refused as a PROPERTY, not a method");
        end

        function directCarriesTheBoxItInherits(tc)
            % The other direction: the one subclass that exists satisfies the
            % widened contract, and the box it exposes is a rid.gen.Constraints.
            gen = GeneratorContractTest.generator(5);
            tc.verifyClass(gen.constraints, "rid.gen.Constraints");
        end

        function theBaseClassCannotBeInstantiated(tc)
            % An abstract class that can be constructed is not abstract, and a
            % Setup argument typed against it would then accept a placeholder that
            % implements nothing.
            tc.verifyError(@() rid.gen.Generator(), ...
                "MATLAB:class:abstract");
        end

        function anIncompleteImplementationCannotBeInstantiated(tc)
            % THE HALF THAT MATTERS. A base class whose abstract block is empty
            % would pass every test above; only refusing an incomplete subclass
            % shows that the contract is enforced by the language rather than
            % described in a comment.
            %
            % Fails if: the abstract block is emptied, or `applied` -- the applied
            % unit the loop commits without assuming its height -- is
            % declared concrete on the base class with a "first sample" default,
            % which would silently make every parametrization a step generator.
            tc.verifyError(@() IncompleteGeneratorDouble(), ...
                "MATLAB:class:abstract");
        end

        function theAppliedUnitIsTheParametrizationsToDeclareNotTheLoops(tc)
            % +rhc commits exactly what `applied` returns. For Direct that
            % is one sample out of a horizon of five, and the remaining four are
            % look-ahead that is applied to nothing.
            %
            % Fails if: `applied` returns the whole trajectory, or a fixed number
            % of samples unrelated to the parametrization's declared unit.
            gen  = GeneratorContractTest.generator(5);
            uSeq = (1:5).' * 0.1;
            uApp = gen.applied(uSeq);
            tc.verifyEqual(size(uApp), [1, gen.nInputs()]);
            tc.verifyEqual(uApp, uSeq(1, :));
        end

    end

    methods (Static, Access = private)

        function gen = generator(nHorizon)
            %GENERATOR  The default SISO parametrization over a small ROI.
            roi = struct("u1_k", [-1 1], "du1_k", [-0.5 0.5], "y_km1", [-2 2]);
            rep = rid.rep.Narx(2, 1);
            gen = rid.gen.Direct(rid.gen.Constraints(rep, roi), nHorizon);
        end

    end

end
