#!/usr/bin/env python3
"""Generate golden test vectors for the rhcInputDesigner criterion.

INDEPENDENT REFERENCE IMPLEMENTATION.

Written, 3 and 5 -- NOT translated from the MATLAB code in
`src/+rid/+crit`, which did not exist when this was written. That independence is the whole point:
two implementations of the same specification can only agree by both being right. Nothing here may
be "corrected" by consulting the MATLAB implementation; if the two disagree, one of them is wrong
and the disagreement is the finding.

Pure standard library on purpose -- no numpy -- so that no library semantics are shared with the
implementation under test.

Usage:  python3 generate_golden.py            (writes *.json beside this script)
"""

import json
import math
import hashlib
import random
import datetime
import pathlib

SPEC = "the method description"


# ---------------------------------------------------------------- the criterion


def refinement(relevance, complexity, noiseLevel, beta, betaNoise):
    """the method description: phi_j = min((a_j c_j)^-beta * v_j^-betaNoise, 1)."""
    out = []
    for a, c, v in zip(relevance, complexity, noiseLevel):
        if a <= 0 or c <= 0 or v <= 0:
            raise ValueError("fields must be strictly positive")
        out.append(min((a * c) ** (-beta) * v ** (-betaNoise), 1.0))
    return out


def demanded_resolution(fields, epsReq0, beta, betaNoise, mode, fillDistanceQ):
    """the method description: epsReq = epsReq0 * phi, floored in certification mode."""
    phi = refinement(fields["relevance"], fields["complexity"], fields["noiseLevel"], beta, betaNoise)
    epsReq = [epsReq0 * p for p in phi]
    if mode == "certification":
        if epsReq0 < fillDistanceQ:
            raise ValueError("rid:crit:uncertifiableBaseline")
        epsReq = [max(e, fillDistanceQ) for e in epsReq]
    return phi, epsReq


def nearest_distances(Znorm, Q):
    """Raw, unnormalized nearest-neighbour distance per supporting point."""
    out = []
    for q in Q:
        if not Znorm:
            out.append(math.inf)
            continue
        best = math.inf
        for z in Znorm:
            # Multiplication, not `** 2`: IEEE-754 requires +, * and sqrt to be
            # correctly rounded, but NOT pow. `(zi - qi) ** 2` goes through libm
            # pow, whose last-bit behaviour is platform-dependent, and that is what
            # made random3d.json irreproducible across machines. With a fixed
            # summation order this expression is now bit-identical on every
            # conforming platform. Numerically inert where pow is already exact.
            d = math.sqrt(sum((zi - qi) * (zi - qi) for zi, qi in zip(z, q)))
            if d < best:
                best = d
        out.append(best)
    return out


def criterion(Znorm, Q, epsReq, psiMax, pNorm):
    """the method description: dTilde, the clipped penalty psi, and J."""
    dRaw = nearest_distances(Znorm, Q)
    dTilde = [d / e for d, e in zip(dRaw, epsReq)]
    psi = [min(max(dt - 1.0, 0.0), psiMax) for dt in dTilde]
    J = (sum(p ** pNorm for p in psi) / len(psi)) ** (1.0 / pNorm)
    return dRaw, dTilde, psi, J


# ---------------------------------------------------------------- case building

CASES = []


def case(name, description, Znorm, Q, epsReq0, *, relevance=None, complexity=None, noiseLevel=None,
         beta=0.5, betaNoise=None, psiMax=1.0, pNorm=1.0, mode="budget", fillDistanceQ=0.0,
         omit_arrays=False):
    nSup = len(Q)
    ones = [1.0] * nSup
    fields = {
        "relevance": relevance if relevance is not None else ones,
        "complexity": complexity if complexity is not None else ones,
        "noiseLevel": noiseLevel if noiseLevel is not None else ones,
    }
    bn = beta if betaNoise is None else betaNoise
    phi, epsReq = demanded_resolution(fields, epsReq0, beta, bn, mode, fillDistanceQ)
    dRaw, dTilde, psi, J = criterion(Znorm, Q, epsReq, psiMax, pNorm)

    expected = {"refinement": phi, "epsReq": epsReq, "J": J}
    if not omit_arrays:
        expected.update({"dRaw": dRaw, "dTilde": dTilde, "psi": psi})

    CASES.append({
        "name": name,
        "description": description,
        "inputs": {
            "Znorm": Znorm, "Q": Q,
            "relevance": fields["relevance"], "complexity": fields["complexity"],
            "noiseLevel": fields["noiseLevel"],
            "epsReq0": epsReq0, "beta": beta, "betaNoise": bn,
            "psiMax": "Inf" if math.isinf(psiMax) else psiMax,
            "pNorm": pNorm, "mode": mode, "fillDistanceQ": fillDistanceQ,
        },
        "expected": expected,
    })


def grid1d(n):
    return [[i / (n - 1)] for i in range(n)]


def grid2d_centred(n):
    h = 1.0 / n
    return [[(i + 0.5) * h, (j + 0.5) * h] for i in range(n) for j in range(n)]


# 1 -- hand-checkable: one sample, one supporting point
case("singleSampleSingleSupport",
     "One sample, one supporting point. dRaw = 0.3, epsReq = 0.1, dTilde = 3, psi = min(2, 1) = 1.",
     [[0.2]], [[0.5]], 0.1)

# 2 -- hand-checkable: 1-D uniform grid, one sample at the centre
case("grid1dOneSample",
     "Eleven supporting points on a 1-D uniform grid, one sample at 0.5, epsReq = 0.25.",
     [[0.5]], grid1d(11), 0.25)

# 3 -- every supporting point resolved: J must be exactly zero (P4)
case("allResolved",
     "Samples placed exactly on the supporting points: every dTilde = 0, so J = 0 exactly.",
     grid1d(5), grid1d(5), 0.25)

# 4 -- everything clipped: J = psiMax (P5)
case("allClipped",
     "A single distant sample: every supporting point saturates, so J = psiMax.",
     [[10.0, 10.0]], grid2d_centred(4), 0.05)

# 5 -- mixed regimes in one case
case("mixedRegimes",
     "Samples chosen so that some supporting points are resolved, some graded, some clipped.",
     [[0.125, 0.125], [0.375, 0.375]], grid2d_centred(4), 0.2)

# 6 -- no samples at all
case("emptySamples",
     "No samples: every dRaw is infinite, so every psi saturates and J = psiMax. Infinite "
     "intermediates are omitted rather than encoded in JSON.",
     [], grid1d(5), 0.2, omit_arrays=True)

# 7 -- fields above and below nominal; the clip at 1 must bite for the below-nominal ones
case("fieldsNonNominal",
     "Fields above nominal refine; fields below nominal are clipped at refinement = 1.",
     [[0.3], [0.7]], grid1d(5), 0.2,
     relevance=[1.0, 2.0, 1.0, 0.5, 1.0],
     complexity=[1.0, 1.0, 4.0, 1.0, 0.25],
     noiseLevel=[1.0, 1.0, 1.0, 1.0, 3.0])

# 8 -- betaNoise differs from beta
case("betaNoiseDiffers",
     "betaNoise != beta: the noise field enters with its own exponent.",
     [[0.3], [0.7]], grid1d(5), 0.2,
     noiseLevel=[1.0, 2.0, 4.0, 1.0, 8.0], beta=0.5, betaNoise=0.25)

# 9 -- pNorm variations
case("pNormTwo", "pNorm = 2 instead of the default 1.",
     [[0.125, 0.125], [0.375, 0.375]], grid2d_centred(4), 0.2, pNorm=2.0)
case("pNormLarge", "pNorm = 32, approaching the worst-case limit.",
     [[0.125, 0.125], [0.375, 0.375]], grid2d_centred(4), 0.2, pNorm=32.0)

# 10 -- unclipped, distance-based limit
case("psiMaxInf",
     "psiMax = Inf: no clipping, so J = mean(dRaw ./ epsReq) - 1, the distance-based limit "
     ".",
     [[0.125, 0.125], [0.375, 0.375]], grid2d_centred(4), 0.2, psiMax=math.inf)

# 11 -- certification floor, inactive then active
case("floorInactive",
     "Certification mode with fillDistanceQ below every refined resolution: the floor never acts.",
     [[0.3], [0.7]], grid1d(5), 0.2,
     complexity=[1.0, 2.0, 4.0, 1.0, 1.0], mode="certification", fillDistanceQ=0.02)
case("floorActive",
     "Certification mode with fillDistanceQ above the most refined resolutions: those are projected "
     "onto the floor, while the baseline itself stays certifiable.",
     [[0.3], [0.7]], grid1d(5), 0.2,
     complexity=[1.0, 2.0, 4.0, 1.0, 1.0], mode="certification", fillDistanceQ=0.12)

# 12 -- seeded pseudo-random cases exercising the vectorized path
rng = random.Random(20260828)
for nDim in (2, 3):
    Q = [[rng.random() for _ in range(nDim)] for _ in range(50)]
    Z = [[rng.random() for _ in range(nDim)] for _ in range(20)]
    comp = [1.0 + 3.0 * rng.random() for _ in range(50)]
    noise = [0.5 + 2.0 * rng.random() for _ in range(50)]
    case("random%dd" % nDim,
         "Fifty supporting points and twenty samples in %d dimensions, with non-nominal complexity "
         "and noise fields. Coordinates are stored explicitly, so the case does not depend on any "
         "random number generator." % nDim,
         Z, Q, 0.15, complexity=comp, noiseLevel=noise, beta=0.5, betaNoise=0.5)


# ---------------------------------------------------------------- writing

here = pathlib.Path(__file__).resolve().parent
digest = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:16]
stamp = datetime.date.today().isoformat()

index = []
for c in CASES:
    c["source"] = {
        "generatedFrom": SPEC,
        "generator": "generate_golden.py",
        "generatorSha256Prefix": digest,
        "date": stamp,
        "note": "Independent Python reference. NOT derived from the MATLAB implementation.",
    }
    path = here / ("%s.json" % c["name"])
    path.write_text(json.dumps(c, indent=2, allow_nan=False) + "\n")
    index.append((c["name"], c["expected"]["J"], len(c["inputs"]["Q"]), len(c["inputs"]["Znorm"])))

print("%-28s %-22s %6s %7s" % ("case", "J", "nSup", "nSamp"))
for name, J, nsup, nsamp in index:
    print("%-28s %-22.17g %6d %7d" % (name, J, nsup, nsamp))
print("\n%d cases written to %s" % (len(CASES), here))
