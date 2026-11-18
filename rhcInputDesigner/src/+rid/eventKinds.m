function kinds = eventKinds()
%EVENTKINDS  The closed set of run-record event kinds, in section 1's order.
%
%       kinds = rid.eventKinds()
%
%   ONE SET, ONE OWNER. `kind` is drawn from a CLOSED set and code
%   branches on it, never on a parsed string; a closed set with two copies is a
%   closed set that can drift, and the drift is invisible until a record written
%   by one is read by the other.
%
%   WHY IT LIVES AT THE TOP LEVEL and not inside a package. It is read by
%   `+rhc`, which emits the events and closes the set by refusing anything
%   outside it, and by `+viz`, which draws them. Putting it in either would make
%   one depend on the other, and `+viz` depending on `+rhc` is the dependency
%   the architecture forbids outright (`the project's coding standard: `+viz` depends on results
%   only). `rid.defaultOptions` and `rid.assertHistory` already establish the
%   top level as where something both sides need belongs -- and that this list
%   had to be written twice is precisely because it had no such home.
%
%   MEMBERSHIP IS NOT A WHITELIST FOR CONSUMERS. `+rhc` uses the set as a guard,
%   so nothing outside it is ever logged; `+viz` uses it only to assign a colour
%   per kind, and draws an unrecognized kind in grey rather than dropping it. A
%   kind added to the specification and emitted before this list catches up is
%   therefore visible on the plot, which is the behaviour that matters at the
%   moment somebody needs to see it.
%
%   THE ORDER IS LOAD-BEARING: consumers colour by index into this array, so
%   reordering repaints every kind after the change.
%
%   Outputs:
%       kinds  (1 x 7) string
kinds = ["surrogateUpdate", "resolutionChange", "coverageReset", ...
         "supportChange", "normalizationChange", "modeChange", "seed"];

end
