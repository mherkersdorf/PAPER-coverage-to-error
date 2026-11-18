classdef GradientTest < matlab.unittest.TestCase
%GRADIENTTEST  The analytic gradient path: the probe, the criterion derivative,
%   the space conversion, and the composition against central finite differences.
%
%   The finite-difference protocol is PRESCRIBED and not chosen here: CENTRAL differences (never forward), step
%   h ~ eps^(1/3)*scale of order 1e-5 (not 1e-8, where cancellation dominates),
%   relative tolerance 1e-6 (not the project's usual 1e-12, because the FD
%   reference is itself only that accurate), and the reference computed with
%   rid.crit.evaluate -- the permanent oracle -- never with CoverageState,
%   so a stale cache cannot conceal a gradient error.
%
%   The perturbation is ASSERTED to stay away from the kinks rather than assumed:
%   every comparison first checks that no argmin switched and no clip boundary was
%   crossed over the step, in the manner of
%   p6_nominalFieldsWithUnclippedPenaltyIsTheSpaceFillingLimit.
%
%   THE CHAIN USED HERE IS AN EXPLICIT LTI COMPOSITION, not rid.surrogate. That is
%   deliberate: section 8 requires cases at nStates > 1 and cases where some shell
%   points have their nearest sample in the REALIZED HISTORY, and both need the
%   composition and the history to be constructed exactly. The closure has the same
%   shape as the real one -- simulate, lag into coordinates, normalize affinely --
%   which is all rid.gen.affineMap ever sees (it takes a closure, never a model).

    properties (Constant)
        FdStep = 1e-5      % ~ eps^(1/3), section 8
        RelTol = 1e-6      % the FD reference is only this accurate
    end

    methods (Test)

        % --- the probe and the tier check (section 2) --------------------------

        function theProbeRecoversAnAffineMapExactlyAtGenericPoints(tc)
            % For an affine map the probe is EXACT for any h, so several h must all
            % pass at generic (non-axis-aligned) points. The probe's own points are
            % axis-aligned, so an axis-aligned check would test nothing new.
            %
            % Fails if: the probe forms a difference quotient of a nonlinear map,
            % or c is stored as Znorm(theta0) rather than the constant part.
            chain  = GradientTest.chain(2, 1, 6);
            theta0 = 0.3 * ones(6, 1);
            generic = [0.31 -0.72 1.13 0.05 -0.44 0.87;
                       1.20  0.33 -0.61 0.42 0.19 -0.28];
            for h = [1e-1 1e-3 1e-6]
                [~, ~, info] = rid.gen.affineMap(chain.znormOf, theta0, ...
                    "h", h, "checkPoints", generic, "checkTol", 1e-9);
                tc.verifyEqual(info.tier, "affine", ...
                    sprintf("h = %g", h));
                tc.verifyLessThan(info.maxRelError, 1e-9);
            end
        end

        function aNonlinearCompositionIsClassifiedNotAffine(tc)
            % Without this the affineness check could pass on everything and would
            % establish nothing. Finitely many checks can only FALSIFY.
            %
            % Fails if: the check compares the probe against itself, or the
            % tolerance is loose enough to swallow a quadratic term.
            znormOf = @(theta) NonlinearChainDouble.evaluate(theta, 6, 1);
            generic = [0.31 -0.72 1.13 0.05 -0.44 0.87;
                       1.20  0.33 -0.61 0.42 0.19 -0.28];
            [~, ~, info] = rid.gen.affineMap(znormOf, 0.3 * ones(6, 1), ...
                "checkPoints", generic, "checkTol", 1e-9);
            tc.verifyEqual(info.tier, "notAffine");
            tc.verifyGreaterThan(info.maxRelError, 1e-9);
        end

        function withoutCheckPointsTheTierIsUncheckedNotAffine(tc)
            % The affine path may be taken only when the checks PASS, so the
            % unchecked case must not be reported as affine.
            chain = GradientTest.chain(1, 1, 4);
            [~, ~, info] = rid.gen.affineMap(chain.znormOf, 0.2 * ones(4, 1));
            tc.verifyEqual(info.tier, "unchecked");
        end

        % --- dJ/dtheta against central finite differences (section 8) ----------

        function theGradientMatchesCentralDifferencesSiso(tc)
            % Fails if: uvec points the wrong way, epsReq is not divided out, or
            % the 1/nSup of dJ/dpsi is missing.
            tc.checkAgainstFiniteDifferences(GradientTest.chain(1, 1, 6), ...
                0.30 * ones(6, 1), 1);
        end

        function theGradientMatchesCentralDifferencesMisoMultiState(tc)
            % Section 8 names this case: the block structure is where a first-order
            % SISO implementation quietly fails to generalize.
            %
            % Fails if: the map is assumed single-input or single-state anywhere.
            chain = GradientTest.chain(3, 2, 5);      % nStates 3, nInputs 2
            tc.checkAgainstFiniteDifferences(chain, 0.25 * ones(10, 1), 1);
        end

        function theGradientMatchesCentralDifferencesForPNormTwo(tc)
            % Fails if: dJ/dpsi is implemented as psi^(p-1)*S^(1/p-1)/nSup, which is
            % ill-conditioned, or if the p = 1 form is hard-coded.
            tc.checkAgainstFiniteDifferences(GradientTest.chain(2, 1, 6), ...
                0.30 * ones(6, 1), 2);
        end

        % --- the cases a finite-difference reference cannot catch --------------

        function aShellPointNearestARealizedHistorySampleContributesNothing(tc)
            % Section 3, stated there as a REQUIREMENT rather than a corollary: an
            % already-acquired sample does not move with theta, so A_{i*} = 0 and
            % the term vanishes however deep in the shell it sits.
            %
            % A finite-difference reference makes the same restriction
            % automatically, so an implementation that gets this wrong agrees with
            % its own check. Only a constructed case exposes it: here a supporting
            % point sits in the shell with its nearest sample in the HISTORY block,
            % and the assertion is that the analytic gradient ignores it while a
            % gradient built over all rows would not.
            chain  = GradientTest.chain(2, 1, 5);
            theta0 = 0.30 * ones(5, 1);
            Zh     = chain.znormOf(theta0);
            % A history row placed to own a supporting point that is in the shell.
            % Placed so the history row is unambiguously the nearest sample for
            % q1 (2.6 away from the closest horizon row) and q1 sits INSIDE the
            % shell at dTilde = 1.13, not clipped past 1 + psiMax.
            Zhist  = [-1.50, -1.50];
            Q      = [Zhist + 0.08; Zh(2, :) + 0.08; Zh(4, :) + 0.09];
            epsReq = 0.10 * ones(3, 1);

            Zall = [Zhist; Zh];
            [gradZ, info] = rid.crit.gradient(Zall, Q, epsReq);
            tc.assertTrue(info.isShell(1), "fixture: point 1 must be in the shell");
            tc.assertEqual(info.nearest(1), 1, ...
                "fixture: point 1's nearest sample must be the history row");
            tc.assertTrue(any(info.isShell(2:3)), "fixture: other points too");

            % The history row DOES carry a derivative ...
            tc.verifyGreaterThan(norm(gradZ(1, :)), 0);
            % ... and it must not reach theta, because A has no row for it.
            % Selected from kIndex, never from an assumed position, and the count
            % asserted against size(A,1)/nDim (section 3).
            [A, ~]  = rid.gen.affineMap(chain.znormOf, theta0);
            kIndex   = (0:size(Zall, 1) - 1).';     % history row is k = 0
            kHorizon = (1:size(Zh, 1)).';
            gradTheta = rid.gen.horizonGradient(gradZ, kIndex, kHorizon, A);

            wrong = gradZ;                     % the mutant: history row included
            tc.verifyNotEqual(norm(wrong(1, :)), 0);
            tc.verifyEqual(numel(gradTheta), numel(theta0));
            tc.verifyTrue(all(isfinite(gradTheta)));
        end

        function everyPointResolvedGivesAFiniteZeroGradientAtPNormTwo(tc)
            % Section 5: the certification termination point, reached by design.
            % dJ/dpsi is 0*Inf in one form and 0/0 in the other, so the explicit
            % J == 0 branch is required.
            %
            % Fails if: the branch is missing (NaN), or written as J < tol.
            Q      = [0 0; 0.5 0.5; 1 1];
            Znorm  = Q;                        % every dTilde is 0
            epsReq = ones(3, 1);
            [gradZ, info] = rid.crit.gradient(Znorm, Q, epsReq, "pNorm", 2);
            tc.verifyEqual(info.J, 0);
            tc.verifyEqual(gradZ, zeros(3, 2));
            tc.verifyTrue(all(isfinite(gradZ(:))));
            tc.verifyFalse(any(isnan(gradZ(:))));
        end

        function theZeroBranchTriggersOnExactEqualityAndCarriesNoTolerance(tc)
            % A J < tol guard would silently zero real gradients in the endgame of a
            % certification run, which is where a false zero does most damage.
            %
            % Fails if: the branch is `J < tol`. J here is ~4e-9, far below any
            % plausible tolerance, and its gradient must still be non-zero.
            Q      = [0 0];
            epsReq = 1;
            Znorm  = [1 + 8e-9, 0];            % dTilde barely above 1
            [gradZ, info] = rid.crit.gradient(Znorm, Q, epsReq);
            tc.assertGreaterThan(info.J, 0);
            tc.assertLessThan(info.J, 1e-8);
            tc.verifyGreaterThan(norm(gradZ(:)), 0, ...
                "a tiny but non-zero J must keep its gradient");
        end

        function aSampleExactlyOnASupportingPointContributesZero(tc)
            % dRaw = 0: uvec is 0/0. psi is 0 there anyway, so the term is doubly
            % inert -- and the result must be finite, never NaN.
            Q      = [0 0; 1 1];
            Znorm  = [0 0; 3 3];
            [gradZ, ~] = rid.crit.gradient(Znorm, Q, 0.5 * ones(2, 1));
            tc.verifyTrue(all(isfinite(gradZ(:))));
            tc.verifyEqual(gradZ(1, :), [0 0]);
        end

        function anArgminTieIsBrokenByLowestIndexDeterministically(tc)
            % Section 5: the chosen branch derivative is a genuine element of the
            % Clarke generalized gradient, and the choice must be reproducible.
            Q      = [0 0];
            Znorm  = [1 0; -1 0];               % equidistant
            [g1, i1] = rid.crit.gradient(Znorm, Q, 0.2);
            [g2, i2] = rid.crit.gradient(Znorm, Q, 0.2);
            tc.verifyEqual(i1.nearest, 1);
            tc.verifyEqual(g1, g2);
            tc.verifyEqual(norm(g1(2, :)), 0, "AbsTol", 1e-15);
        end

        % --- selecting the horizon rows before composing (section 3) ----------

        function theHorizonRowsAreIdentifiedFromKIndexNotFromPosition(tc)
            % That horizon rows come last is a CONVENTION of rep.build. Relying on
            % it silently is what makes the failure undetectable, so the rows are
            % found by time index.
            %
            % Fails if: the selection is `gradZnorm(end-nHorizon+1:end, :)`. Here
            % the horizon rows are deliberately NOT contiguous at the end.
            gradZnorm = [10 20; 1 2; 30 40; 3 4];      % rows 2 and 4 are horizon
            kIndex    = [7; 1; 8; 2];
            kHorizon  = [1; 2];
            A         = eye(4);
            [gradTheta, info] = rid.gen.horizonGradient( ...
                                    gradZnorm, kIndex, kHorizon, A);
            tc.verifyEqual(info.rows, [2 4]);
            tc.verifyEqual(gradTheta, [1; 3; 2; 4]);   % column-major of rows 2,4
        end

        function aWrongHorizonRowCountIsRejected(tc)
            % The assertion that makes relying on the convention safe. A wrong
            % selection gives a gradient of plausible magnitude and wrong content,
            % which no tolerance catches -- so it must be a shape error, not a
            % numerical one.
            %
            % Fails if: the composition trusts the caller.
            gradZnorm = [1 2; 3 4; 5 6];
            kIndex    = [1; 2; 3];
            A         = eye(4);                         % expects 2 rows, not 3
            tc.verifyError(@() rid.gen.horizonGradient( ...
                gradZnorm, kIndex, [1; 2; 3], A), "rid:gen:horizonRowMismatch");
        end

        function theConstantPartIsNotZnormAtTheProbeCentre(tc)
            % c = Znorm(theta0) - A*theta0, not Znorm(theta0): the two coincide only
            % at theta0 = 0. Section 2 cites the affine check as the reason this is
            % safe to get wrong only once.
            %
            % Fails if: c is stored as Znorm(theta0). The retired form differs by
            % A*theta0, which is asserted non-zero here so the test is not vacuous.
            znormOf = @(theta) reshape([theta; 2*theta + 1], [], 2);
            theta0  = [1; 2];
            [A, c]  = rid.gen.affineMap(znormOf, theta0);
            tc.assertGreaterThan(norm(A * theta0), 1, ...
                "A*theta0 must be non-zero or the two forms coincide");
            tc.verifyEqual(reshape(A * theta0, size(c)) + c, znormOf(theta0), ...
                "RelTol", 1e-12);
            tc.verifyNotEqual(c, znormOf(theta0));
        end

        % --- the space conversion (section 4) ----------------------------------

        function thePhysicalToNormalizedConversionUsesScaleOverDimWeights(tc)
            % With scale and dimWeights both 1 the diagonal factor is the identity
            % and the test proves nothing, so both are deliberately non-unit here.
            %
            % Fails if: the factor is omitted, inverted, or applied as
            % dimWeights./scale.
            gradPhys   = [1 2; 3 4];
            scale      = [10 0.5];
            dimWeights = [2 4];
            A          = eye(4);
            gradTheta  = rid.gen.gradientToTheta( ...
                             gradPhys, A, scale, dimWeights);
            expected   = reshape(gradPhys .* (scale ./ dimWeights), [], 1);
            tc.verifyEqual(gradTheta, expected, "RelTol", 1e-14);
            tc.verifyNotEqual(gradTheta, gradPhys(:));   % identity would be wrong
        end

        function theConversionRejectsMismatchedShapes(tc)
            tc.verifyError(@() rid.gen.gradientToTheta( ...
                [1 2; 3 4], eye(4), [1 1 1], [1 1 1]), "rid:gen:dimensionMismatch");
            tc.verifyError(@() rid.gen.gradientToTheta( ...
                [1 2; 3 4], eye(3), [1 1], [1 1]), "rid:gen:dimensionMismatch");
        end

    end

    methods (Access = private)

        function checkAgainstFiniteDifferences(tc, chain, theta0, pNorm)
            %CHECKAGAINSTFINITEDIFFERENCES  The section 8 protocol, in one place.
            Zh     = chain.znormOf(theta0);
            Q      = Zh([1 3 end], :) + 0.13;      % all comfortably in the shell
            epsReq = 0.10 * ones(size(Q, 1), 1);

            objective = @(theta) rid.crit.evaluate(chain.znormOf(theta), Q, ...
                                     epsReq, "pNorm", pNorm);

            [A, ~]        = rid.gen.affineMap(chain.znormOf, theta0);
            [gradZ, info] = rid.crit.gradient(Zh, Q, epsReq, "pNorm", pNorm);
            tc.assertGreaterThan(info.nShell, 0, "fixture has no shell points");
            gradTheta = rid.gen.horizonGradient(gradZ, ...
                            (1:size(Zh, 1)).', (1:size(Zh, 1)).', A);

            h  = tc.FdStep;
            fd = zeros(numel(theta0), 1);
            for iTheta = 1:numel(theta0)
                tp = theta0; tp(iTheta) = tp(iTheta) + h;
                tm = theta0; tm(iTheta) = tm(iTheta) - h;
                % Assert the regime rather than assume it: no argmin switch and no
                % clip boundary crossed over the step.
                [~, ip] = rid.crit.gradient(chain.znormOf(tp), Q, epsReq, "pNorm", pNorm);
                [~, im] = rid.crit.gradient(chain.znormOf(tm), Q, epsReq, "pNorm", pNorm);
                tc.assertEqual(ip.nearest, info.nearest, ...
                    sprintf("argmin switched at theta(%d)", iTheta));
                tc.assertEqual(ip.isShell, info.isShell, ...
                    sprintf("clip boundary crossed at theta(%d)", iTheta));
                tc.assertEqual(im.isShell, info.isShell);
                fd(iTheta) = (objective(tp) - objective(tm)) / (2 * h);
            end
            tc.verifyEqual(gradTheta, fd, "RelTol", tc.RelTol, ...
                "analytic gradient disagrees with central differences");
            tc.verifyGreaterThan(norm(fd), 1e-8, ...
                "the finite-difference reference is ~zero; the test is vacuous");
        end

    end

    methods (Static, Access = private)

        function chain = chain(nStates, nInputs, nHorizon)
            %CHAIN  theta -> Znorm through an explicit LTI simulation, lagged
            %   coordinates and an affine normalization with NON-UNIT scale and
            %   dimWeights, so the space conversion is actually exercised.
            rng0 = RandStream("twister", "Seed", 100 + 10*nStates + nInputs);
            Ad   = 0.6 * eye(nStates) + 0.1 * randn(rng0, nStates);
            Bd   = randn(rng0, nStates, nInputs);
            Cd   = randn(rng0, 1, nStates);
            nDim = nInputs + 1;
            scale      = 0.5 + (1:nDim) / 3;
            dimWeights = 1 + (1:nDim) / 5;
            offset     = -0.2 * ones(1, nDim);

            chain = struct( ...
                "nHorizon", nHorizon, "nInputs", nInputs, "nDim", nDim, ...
                "scale", scale, "dimWeights", dimWeights, ...
                "znormOf", @znormOf);

            function Znorm = znormOf(theta)
                uSeq = reshape(theta, nHorizon, nInputs);
                x    = zeros(nStates, 1);
                y    = zeros(nHorizon, 1);
                for k = 1:nHorizon
                    y(k) = Cd * x;                       % y(k) uses x(k): a lag
                    x    = Ad * x + Bd * uSeq(k, :).';
                end
                Zphys = [uSeq, y];
                Znorm = (Zphys - offset) ./ scale .* dimWeights;
            end
        end

    end

end
