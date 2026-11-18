classdef NominalTest < matlab.unittest.TestCase
%NOMINALTEST  The v1 field estimator, and the "one definition of nominal" rule.
%
%   Two things are worth asserting about an estimator that returns ones: that it
%   really is substitutable for a real one, and that its value comes from the
%   single definition rather than from a second literal.

    methods (Test)

        function itReturnsTheNominalValueForEverySupportingPoint(tc)
            % Fails if: the shape is (1 x nSup), or nSup is taken from something
            % other than Q.
            est = rid.fields.Nominal();
            for nSup = [1 7 40]
                Q = rand(nSup, 2);
                [value, info] = est.estimate(NominalTest.hist(), ...
                    rid.surrogate.LinearFirstOrder(0.5, 0.1, 1), ...
                    rid.rep.Narx(2, 1), NominalTest.nrm(), Q);
                tc.verifyEqual(value, ones(nSup, 1));
                tc.verifyEqual(info.nSup, nSup);
                tc.verifyTrue(all(value > 0) && all(isfinite(value)));
            end
        end

        function itAgreesWithTheSingleDefinitionOfNominal(tc)
            % nominalFields is the SOLE constructor
            % of the nominal case, so this estimator must not carry a second
            % literal that could drift from it.
            %
            % Fails if: Nominal writes ones(nSup, 1) of its own. Change the
            % definition and this test catches the copy that did not follow.
            est   = rid.fields.Nominal();
            nSup  = 9;
            value = est.estimate(NominalTest.hist(), ...
                rid.surrogate.LinearFirstOrder(0.5, 0.1, 1), ...
                rid.rep.Narx(2, 1), NominalTest.nrm(), rand(nSup, 2));
            nominal = rid.fields.nominalFields(nSup);
            tc.verifyEqual(value, nominal.relevance);
            tc.verifyEqual(value, nominal.complexity);
            tc.verifyEqual(value, nominal.noiseLevel);
        end

        function nominalFieldsReduceTheDesignToPureSpaceFilling(tc)
            % The property that makes this the v1 default: the refinement product
            % is 1, so epsReq is epsReq0 everywhere and J measures uniform
            % coverage (section 3).
            %
            % Fails if: the estimator returns anything but the nominal value --
            % epsReq would then vary across Q with no field asking it to.
            est    = rid.fields.Nominal();
            Q      = rand(12, 2);
            ones12 = est.estimate(NominalTest.hist(), ...
                rid.surrogate.LinearFirstOrder(0.5, 0.1, 1), ...
                rid.rep.Narx(2, 1), NominalTest.nrm(), Q);
            fields = struct("relevance", ones12, "complexity", ones12, ...
                            "noiseLevel", ones12);
            [epsReq, info] = rid.crit.demandedResolution(fields, "epsReq0", 0.2);
            tc.verifyEqual(epsReq, 0.2 * ones(12, 1), "RelTol", 1e-12);
            tc.verifyEqual(info.fracRefined, 0);
        end

        function itDeclaresItsCalibrationAndThatItIsNotOptimistic(tc)
            % Section 1 item 4: every estimator states what "1" means, and the
            % statement goes into the run log. Section 5: an estimator with no
            % uncertainty cannot be optimistic about it, and must say so rather
            % than leave the field unset.
            %
            % Fails if: calibration is empty, which would put a blank line in the
            % run record where the meaning of the numbers should be.
            est = rid.fields.Nominal();
            tc.verifyClass(est.calibration, "string");
            tc.verifyGreaterThan(strlength(est.calibration), 0);
            tc.verifyFalse(est.isOptimistic);
            tc.verifySubstring(est.describe(), "Nominal");
        end

        function itIsPureAndIgnoresEveryArgumentItIsGiven(tc)
            % Section 7's contract: estimate is pure and deterministic. This
            % estimator asserts something stronger -- that no data could change its
            % answer -- so two wildly different histories must give the same value.
            %
            % Fails if: a future edit makes it read hist or srg and forgets that
            % nominal means nominal.
            est = rid.fields.Nominal();
            Q   = rand(6, 2);
            a = est.estimate(NominalTest.hist(), ...
                rid.surrogate.LinearFirstOrder(0.5, 0.1, 1), ...
                rid.rep.Narx(2, 1), NominalTest.nrm(), Q);
            b = est.estimate(wellFormedHistory(rand(400, 1), "y", rand(400, 1)), ...
                rid.surrogate.LinearFirstOrder(0.9, 0.5, 4), ...
                rid.rep.Narx(3, 2), NominalTest.nrm(), Q);
            tc.verifyEqual(a, b);
        end

        function anEmptySupportSetIsAnError(tc)
            % nSup comes from Q, so an empty Q is a caller error rather than a
            % zero-length field.
            est = rid.fields.Nominal();
            tc.verifyError(@() est.estimate(NominalTest.hist(), ...
                rid.surrogate.LinearFirstOrder(0.5, 0.1, 1), ...
                rid.rep.Narx(2, 1), NominalTest.nrm(), zeros(0, 2)), ...
                "rid:fields:invalidFieldSize");
        end

    end

    methods (Static, Access = private)

        function h = hist()
            h = wellFormedHistory((1:8).', "y", (1:8).');
        end

        function n = nrm()
            n = rid.roi.Normalization([0 1; 0 1], ["u1_k" "y_km1"]);
        end

    end

end
