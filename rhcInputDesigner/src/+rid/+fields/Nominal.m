classdef Nominal
%NOMINAL  The field estimator that asks for no refinement anywhere.
%
%   The v1 default, and the one that reduces the whole design to pure
%   space-filling: with every field at its nominal value the refinement product is
%   1, `epsReq = epsReq0` everywhere, and `J` measures uniform coverage of the ROI.
%
%   IT IS A FIELD ESTIMATOR LIKE ANY OTHER, not a special case the loop knows
%   about. That is the point: the run holds three estimators, one per field, and
%   swapping a real complexity estimator in changes nothing else. An
%   estimator estimates ONE field -- section 1's calibration contract is per field,
%   and section 7's `estimate` returns `(nSup x 1)` -- so the same object stands in
%   for all three until there is something better to put in any of them.
%
%   THE VALUE COMES FROM `rid.fields.nominalFields`, not from a literal here. That
%   function is the sole constructor of the nominal case (section 7), and its three
%   members are the same number by definition, so reading one of them reads the
%   definition. Writing `ones(nSup, 1)` in this file would create the second
%   definition the rule exists to prevent.
%
%   IT IS NOT OPTIMISTIC, and says so. Section 5's trap is that a region with few
%   samples looks simple and quiet whether or not it is, so the loop confirms its
%   own neglect; an optimistic variant lets the uncertainty itself drive refinement
%   into unexplored regions. This estimator has no uncertainty to be optimistic
%   about -- it is not an estimate at all -- so `isOptimistic` is false, and that
%   is a statement about the estimator rather than a defect in it.
%
%   VALUE SEMANTICS and no state: `estimate` ignores every argument, which is what
%   makes it trivially pure and deterministic (section 7's contract).
%
%   Construction:
%       est = rid.fields.Nominal()

    properties (SetAccess = immutable)
        % What "1" means for this estimator, written into the run log (section 1,
        % item 4). Every estimator declares this; an estimator producing only
        % relative values states its convention here rather than hiding it.
        calibration  (1,1) string
        isOptimistic (1,1) logical
    end

    methods

        function obj = Nominal()
            obj.calibration = "nominal by construction: every supporting point " + ...
                "asks for the baseline resolution epsReq0 and no refinement, " + ...
                "so the design reduces to pure space-filling.";
            obj.isOptimistic = false;
        end

        function [value, info] = estimate(obj, hist, srg, rep, nrm, Q) %#ok<INUSD>
            %ESTIMATE  The nominal field: one value per supporting point, all 1.
            %
            %   Every argument is ignored, and deliberately so -- this estimator
            %   asserts that no data could change its answer. The full signature is
            %   kept because it is the interface a real estimator implements
            %   (section 7), and a shorter one here would make Nominal
            %   un-substitutable.
            %
            %   Inputs:
            %       hist  (1x1) struct  the signal history
            %       srg   (1x1) object  the current surrogate
            %       rep   (1x1) object  the representation
            %       nrm   (1x1) object  the normalization
            %       Q     (nSup x nDim) supporting points, normalized
            %   Outputs:
            %       value (nSup x 1)    all ones; finite and strictly positive
            %       info  (1x1) struct  .calibration .isOptimistic .nSup
            arguments
                obj  (1,1) rid.fields.Nominal
                hist (1,1) struct
                srg  (1,1)
                rep  (1,1)
                nrm  (1,1)
                Q    (:,:) double {mustBeReal, mustBeFinite}
            end
            if size(Q, 1) < 1
                error("rid:fields:invalidFieldSize", ...
                    "Q must have at least one supporting point.");
            end
            % The single definition of nominal; see the class header.
            value = rid.fields.nominalFields(size(Q, 1)).relevance;
            info  = struct("calibration",  obj.calibration, ...
                           "isOptimistic", obj.isOptimistic, ...
                           "nSup",         size(Q, 1));
        end

        function s = describe(obj)
            %DESCRIBE  One line for the run record.
            s = sprintf("rid.fields.Nominal (isOptimistic = %d)", obj.isOptimistic);
        end

    end

end
