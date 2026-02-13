"""
colors.py — A utility for picking visually distinct colors from the ROOT
color wheel, ordered for good contrast when cycling through multiple
histograms.

ROOT colors are referenced by integer codes. The base hues and their
canonical numbers are:
    kRed=632, kGreen=416, kBlue=600, kYellow=400, kMagenta=616,
    kCyan=432, kOrange=800, kTeal=840, kViolet=880, kPink=900,
    kSpring=820, kAzure=860

Each base can be offset by -10..+4 to get lighter/darker or more
saturated variants (e.g. kRed-4, kBlue+2).  Positive offsets go
darker/more saturated; negative offsets go lighter/paler.

The list below is ordered so that consecutive entries are drawn from
opposite sides of the wheel, maximising contrast between adjacent
histograms in a stack or overlay.
"""

from itertools import cycle
from typing import Iterator

# ---------------------------------------------------------------------------
# ROOT color name → integer mappings (from Rtypes.h / TColor.h)
# ---------------------------------------------------------------------------
kRed     = 632
kGreen   = 416
kBlue    = 600
kYellow  = 400
kMagenta = 616
kCyan    = 432
kOrange  = 800
kTeal    = 840
kViolet  = 880
kPink    = 900
kSpring  = 820
kAzure   = 860

# ---------------------------------------------------------------------------
# Core palette — ROOT integer color codes
# ---------------------------------------------------------------------------
# Consecutive pairs are near-complementary (opposite on the wheel) so that
# the 1st and 2nd histogram, 3rd and 4th, etc. always contrast well.
# Offsets are chosen to land on vivid, clean shades rather than the
# sometimes-washed-out bare base values.

COLORS: list[int] = [
    kRed    + 1,   # bright red
    kAzure  + 1,   # rich azure-blue
    kGreen  + 2,   # vivid green
    kOrange + 1,   # warm orange
    kViolet + 1,   # purple-violet
    kYellow - 7,   # golden yellow (base is pale; -7 gives a richer tone)
    kCyan   + 1,   # cyan-teal
    kPink   - 4,   # rose-pink
    kSpring + 4,   # yellow-green
    kMagenta+ 2,   # deep magenta
    kTeal   + 2,   # dark teal
    kRed    - 7,   # soft coral (lighter red variant)
]


# ---------------------------------------------------------------------------
# Access helpers
# ---------------------------------------------------------------------------

def get_color(index: int) -> int:
    """Return a ROOT color code by *index*, wrapping around automatically.

    Usage::

        for i, hist in enumerate(histograms):
            hist.SetLineColor(get_color(i))
            hist.SetFillColor(get_color(i))
    """
    return COLORS[index % len(COLORS)]


def color_cycle() -> Iterator[int]:
    """Return an infinite iterator that cycles through the palette.

    Usage::

        colors = color_cycle()
        for hist in histograms:
            hist.SetLineColor(next(colors))
    """
    return cycle(COLORS)


# ---------------------------------------------------------------------------
# Optional: quick printout of the palette (run this file directly)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("ROOT color palette — integer codes")
    print(f"{'Index':>6}  {'Code':>6}  Label")
    print("-" * 32)
    labels = [
        "kRed+1",    "kAzure+1",  "kGreen+2",  "kOrange+1",
        "kViolet+1", "kYellow-7", "kCyan+1",   "kPink-4",
        "kSpring+4", "kMagenta+2","kTeal+2",   "kRed-7",
    ]
    for i, (code, label) in enumerate(zip(COLORS, labels)):
        print(f"  [{i:2d}]  {code:>6}  {label}")