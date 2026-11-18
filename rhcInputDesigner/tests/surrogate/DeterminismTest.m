classdef DeterminismTest < matlab.unittest.TestCase
%DETERMINISMTEST  The surrogate is deterministic and touches no global RNG state.
%
%   Nothing in +rep or +surrogate is randomized: `simulate` is pure, `build` is a
%   reindexing of the signals, and `update` is least squares. The point of this
%   file is to assert that this stays true.
%
%   The global-stream check is the one that matters. the project's coding standard requires randomness
%   to be confined to entry points, and a library function that quietly consumed
%   global stream state would make every SEEDED run downstream of it irreproducible
%   -- a failure that appears only as an experiment that will not replay, long
%   after the change that caused it.

    properties (Constant)
        TS = 0.05
    end

    methods (Access = private)
        function hist = seededHistory(tc, seed)
            % A local RandStream, never the global one: the test's own randomness
            % must not be the thing under test.
            rs = RandStream("twister", "Seed", seed);
            u  = 2 * rand(rs, 300, 1) - 1;
            srg = rid.surrogate.LinearFirstOrder(0.8, tc.TS, 1.7, "offset", 0.2);
            hist = struct("u", u, "y", srg.simulate(u, srg.steadyState(u(1))), ...
                "tS", tc.TS, "x0", srg.steadyState(u(1)));
        end
    end

    methods (Test)

        function repeatedRunsAreBitIdentical(tc)
            hist = tc.seededHistory(101);
            srg  = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);

            for src = ["measured", "simulated"]
                rep = rid.rep.Narx(2, 1, "pastSource", src);
                [Z1, k1] = rep.build(hist, srg);
                [Z2, k2] = rep.build(hist, srg);
                tc.verifyEqual(Z1, Z2, sprintf("build, pastSource = %s", src));
                tc.verifyEqual(k1, k2);
            end

            f1 = srg.update(hist);
            f2 = srg.update(hist);
            tc.verifyEqual(f1.timeConstant, f2.timeConstant);
            tc.verifyEqual(f1.gain, f2.gain);
            tc.verifyEqual(f1.offset, f2.offset);

            rs = RandStream("twister", "Seed", 102);
            uSeqBatch = rand(rs, 20, 1, 9);
            tc.verifyEqual(srg.simulateBatch(uSeqBatch, 0.3), ...
                           srg.simulateBatch(uSeqBatch, 0.3));
        end

        function noGlobalRandomStreamStateIsConsumed(tc)
            before = rng();

            hist = tc.seededHistory(103);
            srg  = rid.surrogate.LinearFirstOrder(2, tc.TS, 3);
            rep  = rid.rep.Narx(2, 1, "pastSource", "simulated");
            rep.build(hist, srg);
            rep.describe();
            srg.simulate(hist.u, hist.x0);
            srg.simulateBatch(reshape(hist.u, [], 1, 1), hist.x0);
            srg.steadyState(0.5);
            srg.update(hist);
            MinimalSurrogate().steadyState(1);
            sch = rid.surrogate.ImmediateScheduler("updateInterval", 10);
            sch.poll(hist, srg);

            after = rng();
            tc.verifyEqual(after.Type,  before.Type);
            tc.verifyEqual(after.Seed,  before.Seed);
            tc.verifyEqual(after.State, before.State);
        end

        function theSameSeedGivesTheSameFit(tc)
            f1 = rid.surrogate.LinearFirstOrder(2, tc.TS, 3).update(tc.seededHistory(104));
            f2 = rid.surrogate.LinearFirstOrder(2, tc.TS, 3).update(tc.seededHistory(104));
            tc.verifyEqual(f1.describe(), f2.describe());
        end

    end

end
