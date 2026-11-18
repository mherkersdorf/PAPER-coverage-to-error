classdef CriterionContractTest < matlab.unittest.TestCase
%CRITERIONCONTRACTTEST  The defaults of rid.crit.Criterion and CoverageState's overrides.
%
%   The defaults are each the honest answer of a criterion that LACKS a concept,
%   and the loop relies on exactly those answers; CoverageState must override
%   every one of them, or the loop would treat the coverage criterion as one
%   without certification, monotonicity or resolution.

    methods (Test)

        function theDefaultsDescribeACriterionWithoutTheConcepts(tc)
            c = MeanSquareCriterion([0.5 0.5]);

            tc.verifyFalse(c.isSatisfied());
            tc.verifyFalse(c.isMonotone());
            tc.verifyEmpty(c.resolution());
            tc.verifyTrue(isnan(c.nActive()));
            tc.verifyTrue(c.isConsistentWith(rand(3, 2)));
            tc.verifySameHandle(c.withoutQ(), c);
            tc.verifyError(@() c.setResolution(1), "rid:crit:noResolution");
        end

        function coverageStateOverridesEveryDefault(tc)
            Q      = [0 0; 1 0; 0 1; 1 1];
            epsReq = [0.3; 0.3; 0.4; 0.4];
            s      = rid.crit.CoverageState(Q, epsReq, "psiMax", 2);
            s.add([0.1 0.1]);

            tc.verifyClass(s, "rid.crit.CoverageState");
            tc.verifyTrue(isa(s, "rid.crit.Criterion"));
            tc.verifyTrue(s.isMonotone());
            tc.verifyEqual(s.resolution(), epsReq);
            % One point resolved (distance 0.141 <= 0.3), three not.
            tc.verifyEqual(s.nActive(), 3);
            tc.verifyTrue(s.isConsistentWith(Q));
            tc.verifyFalse(s.isConsistentWith(Q + 0.1));
            tc.verifyEqual(size(s.withoutQ().Q), [0 2]);
            s.setResolution([0.1; 0.3; 0.4; 0.4]);
            tc.verifyEqual(s.nActive(), 4, "setResolution must still act");
        end

        function aFreshCopyHasTheConfigurationAndNothingCommitted(tc)
            Q = [0 0; 1 1];
            s = rid.crit.CoverageState(Q, [0.2; 0.5], "psiMax", 3, "pNorm", 2);
            s.add([0.9 0.9]);
            fresh = s.freshCopy();

            tc.verifyNotSameHandle(fresh, s);
            tc.verifyEqual(fresh.Q, Q);
            tc.verifyEqual(fresh.epsReq, [0.2; 0.5]);
            tc.verifyEqual([fresh.psiMax, fresh.pNorm], [3, 2]);
            tc.verifyEqual(fresh.nSamples, 0);
            tc.verifyEqual(fresh.dRaw, [Inf; Inf]);
            tc.verifyEqual(s.nSamples, 1, "the original must be untouched");
        end

        function coverageStateDescribesItself(tc)
            s = rid.crit.CoverageState([0 0; 1 1], [0.2; 0.2], "psiMax", 3.6);

            tc.verifyEqual(s.describe(), "rid.crit.CoverageState: 2 supporting " + ...
                "point(s) in 2 dimension(s), psiMax 3.6, pNorm 1");
        end

    end

end
