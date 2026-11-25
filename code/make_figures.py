#!/usr/bin/env python3
"""Make the figures of the paper: figures/fig1..fig5 as PNG (200 dpi) and SVG.

Data come from the existing check scripts (p3_spindle.py, payoff_green_check.py).
Run from anywhere:  OMP_NUM_THREADS=1 nice -n 19 python3 code/make_figures.py
Every figure goes through check_layout (no text over text, over a line, arrow, marker or
patch outline, or off the canvas); a violation raises and the exit status is nonzero.
Layout check only, no files written (about 2 minutes):  python3 code/make_figures.py --check
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "figures")
sys.path.insert(0, HERE)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Arc

BLUE, ORANGE, TEAL, GREY, INK = "#1f4e79", "#c8553d", "#3b8b8b", "#8a8f98", "#222831"
LCOL = {0: BLUE, 1: ORANGE, 2: TEAL, 3: "#8e6bb0"}
plt.rcParams.update({
    "font.size": 11.5, "axes.labelsize": 12.5, "legend.fontsize": 10.5,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "font.family": "serif", "font.serif": ["CMU Serif", "cmr10", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False, "svg.fonttype": "none",
})


CHECK_ONLY = "--check" in sys.argv   # draw and layout-check only, write no files

# ------------------------------------------------------------ layout check
# Every figure is checked after drawing (see check_layout). A Text is exempt from the
# "text on an artist" test only if it carries a gid starting with "own-label": it is the
# label of a filled box or bar it sits inside on purpose (panel (a) subshell labels and
# the period numbers of panel (b) of the Madelung figure). Exempt texts are still tested
# against every other text and against the canvas edge. Tick labels touching axis spines
# are not tested (spines are not collected).
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle as _Rect
from matplotlib.transforms import Bbox
from matplotlib.legend import Legend


def _tight(bb, pad=0.0):
    return Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)


def check_layout(fig, name):
    """Raise AssertionError if any text overlaps another text, touches a line, arrow,
    marker or patch outline, or leaves the canvas."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi = fig.dpi
    texts = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        if t.axes is not None and t in (t.axes.title, t.axes._left_title, t.axes._right_title) and not t.get_text():
            continue
        bb = t.get_window_extent(r)
        if not np.isfinite([bb.x0, bb.x1, bb.y0, bb.y1]).all() or bb.width == 0:
            continue
        # skip tick labels that matplotlib keeps but does not draw (outside the view)
        texts.append((t, bb))
    # tick labels outside the axes view are not drawn: drop them
    drawn = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for t in axis.get_majorticklabels() + axis.get_minorticklabels():
                drawn.add(id(t))
    inview = []
    for ax in fig.axes:
        ab = ax.get_window_extent(r)
        for axis, hor in ((ax.xaxis, True), (ax.yaxis, False)):
            ticks = list(axis.get_major_ticks()) + list(axis.get_minor_ticks())
            lo, hi = axis.get_view_interval()
            for tk in ticks:
                loc = tk.get_loc()
                if not (min(lo, hi) - 1e-9 <= loc <= max(lo, hi) + 1e-9):
                    tk.label1.set_visible(False)
    texts = [(t, bb) for t, bb in texts if t.get_visible()]
    # artists
    items = []   # (label, path in display coords or None, marker boxes)
    for ax in fig.axes:
        for ln in ax.lines:
            if len(ln.get_xdata()) == 0 or not ln.get_visible():
                continue
            if ln.get_color() == "none":
                continue
            path = ln.get_transform().transform_path(ln.get_path())
            if ln.get_linestyle() in ("None", "none", " ", ""):
                path = None
            boxes = []
            if ln.get_marker() not in ("None", "none", None, ""):
                pts = ln.get_transform().transform(np.column_stack([ln.get_xdata(), ln.get_ydata()]))
                h = ln.get_markersize() * dpi / 72 / 2 + ln.get_markeredgewidth()
                boxes = [Bbox.from_extents(x - h, y - h, x + h, y + h) for x, y in pts if np.isfinite([x, y]).all()]
            items.append((f"line {ln.get_label()!r}", path, boxes))
        for pa in ax.patches:
            if not pa.get_visible() or pa is ax.patch:
                continue
            path = _flat(pa.get_transform().transform_path(pa.get_path()))
            if isinstance(pa, _Rect) and not pa.get_ec()[3] and pa.get_fc()[3] == 0:
                continue
            items.append((f"patch {type(pa).__name__}", path, []))
        for ch in ax.get_children():
            if isinstance(ch, Text) and hasattr(ch, "arrow_patch") and ch.arrow_patch is not None:
                pa = ch.arrow_patch
                items.append(("annotation arrow", _flat(pa.get_path().transformed(pa.get_transform())), []))
        for lg in [ax.get_legend()]:
            if lg is not None and lg.get_frame_on():
                fb = lg.get_frame().get_window_extent(r)
                items.append(("legend frame", __import__("matplotlib").path.Path.unit_rectangle().transformed(
                    __import__("matplotlib.transforms", fromlist=["x"]).BboxTransformTo(fb)), []))
    problems = []
    fb = fig.bbox
    for i, (t, bb) in enumerate(texts):
        s = repr(t.get_text())[:30]
        if bb.x0 < fb.x0 - 0.5 or bb.y0 < fb.y0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y1 > fb.y1 + 0.5:
            problems.append(f"text {s} extends beyond the canvas")
        for t2, bb2 in texts[i + 1:]:
            if bb.overlaps(bb2) and bb.fully_overlaps(bb2) or (bb.overlaps(bb2) and _inter(bb, bb2) > 0.5):
                problems.append(f"text {s} overlaps text {repr(t2.get_text())[:30]}")
        if (t.get_gid() or "").startswith("own-label"):
            continue
        tb = _tight(bb, 0.5)
        for lab, path, boxes in items:
            if path is not None and path.intersects_bbox(tb, filled=False):
                problems.append(f"text {s} touches {lab}")
            elif any(tb.overlaps(b) for b in boxes):
                problems.append(f"text {s} touches a marker of {lab}")
    if problems:
        raise AssertionError(f"layout check FAILED for {name}:\n  " + "\n  ".join(sorted(set(problems))))
    print(f"layout {name}: PASS ({len(texts)} texts, {len(items)} artists)")


def _flat(path):
    """Polyline version of a path (Bezier control points of circles and arcs lie off the curve)."""
    from matplotlib.path import Path
    verts, codes = [], []
    for poly in path.to_polygons(closed_only=False):
        verts.extend(poly)
        codes.extend([Path.MOVETO] + [Path.LINETO] * (len(poly) - 1))
    return Path(np.array(verts), codes) if verts else path


def _inter(a, b):
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return min(w, h)


def save(fig, name):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout(pad=0.6)
    check_layout(fig, name)
    if CHECK_ONLY:
        plt.close(fig)
        return
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".svg"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 1 spindle
def fig1():
    s = np.linspace(0, np.pi, 600)
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    cols = {0.25: TEAL, 0.5: BLUE, 1.0: ORANGE}
    for mu, c in cols.items():
        rho = 2 * mu * np.sin(s)
        lw = 2.6 if mu == 0.5 else 1.8
        ax.plot(s, rho, color=c, lw=lw, label={0.25: r"$\mu=1/4$", 0.5: r"$\mu=1/2$ (round $S^4$)", 1.0: r"$\mu=1$"}[mu])
        ax.plot(s, -rho, color=c, lw=lw)
        # cone tangent at the pole P_0 (radius = 2 mu * distance)
        t = np.linspace(0, 1.1, 2)
        ax.plot(t, 2 * mu * t, color=c, lw=0.9, ls=(0, (3, 3)))
        ax.plot(t, -2 * mu * t, color=c, lw=0.9, ls=(0, (3, 3)))
    ax.axhline(0, color=GREY, lw=0.8)
    ax.set_xlim(-0.3, np.pi + 0.3); ax.set_ylim(-2.5, 2.35)
    ax.set_aspect("equal")
    ax.set_xticks([0, np.pi / 2, np.pi])
    ax.set_xticklabels([r"$0$", r"$\theta=\pi/2$", r"$\pi$"])
    ax.set_xlabel(r"arc length from the nucleus, $\theta$ (units of $\sqrt{R}$)")
    ax.set_ylabel(r"circumferential radius$/\sqrt{R}$" "\n" r"$2\mu\,\sin\theta$")
    ax.text(0.0, -2.38, r"$P_0$", ha="center", color=INK)
    ax.text(np.pi, -2.38, r"$P_\infty$", ha="center", color=INK)
    ax.spines["left"].set_visible(True)
    ax.legend(frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.16), columnspacing=1.2, handlelength=1.6)
    save(fig, "fig1_spindle")


# ---------------------------------------------------------------- 2 S^4 picture
def fig2():
    fig, ax = plt.subplots(figsize=(6.6, 3.5))
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Circle((0, 0), 1, fc="#eef3f8", ec=INK, lw=1.6, zorder=1))
    th = np.linspace(0, np.pi, 200)
    ax.plot(np.cos(th), -0.3 * np.sin(th), color=ORANGE, lw=2.2, zorder=3)          # equator front
    ax.plot(np.cos(th), 0.3 * np.sin(th), color=ORANGE, lw=1.4, ls=(0, (4, 3)), zorder=3)  # back
    # meridian
    ax.plot(0.5 * np.sin(th) * 0 + 0.0 * th, np.linspace(-1, 1, 200), color="none")
    mer = np.linspace(-np.pi / 2, np.pi / 2, 200)
    ax.plot(0.36 * np.cos(mer), np.sin(mer), color=GREY, lw=0.9, zorder=2)
    ax.plot([0, 0], [-1, 1], color=GREY, lw=0.9, ls=(0, (2, 3)), zorder=2)
    for y, lab, dy in ((-1, r"$P_0$  nucleus, $r=0$", -0.12), (1, r"$P_\infty$, $r=\infty$", 0.11)):
        ax.plot(0, y, "o", color=BLUE, ms=7, zorder=5)
        ax.text(0, y + dy, lab, ha="center", va="top" if dy < 0 else "bottom", color=BLUE, fontsize=11.5)
    ax.text(1.04, -0.31, r"inversion sphere $r=R$" "\n" r"($\theta=\pi/2$)", color=ORANGE, va="center", fontsize=10.5)
    ax.text(-0.40, -0.45, r"$r<R$", color=INK, ha="center", fontsize=10)
    ax.text(-0.40, 0.45, r"$r>R$", color=INK, ha="center", fontsize=10)
    # Kelvin swap arrow outside the sphere
    arr = FancyArrowPatch((-1.4, -0.75), (-1.4, 0.75), connectionstyle="arc3,rad=-0.35",
                          arrowstyle="<|-|>", mutation_scale=18, color=TEAL, lw=1.9, zorder=6, shrinkA=6, shrinkB=6)
    arr.set_clip_on(False)
    ax.add_patch(arr)
    ax.text(-2.75, 0.0, "Kelvin inversion\n" r"$r\mapsto R^2/r$" "\n(swaps the poles,\nreflects in the equator)",
            color=TEAL, ha="center", va="center", fontsize=10.5)
    ax.set_xlim(-3.75, 2.3); ax.set_ylim(-1.45, 1.4)
    save(fig, "fig2_s4_picture")


# ---------------------------------------------------------------- 3 Madelung groups
def fig3():
    names = "spdfg"
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.6), gridspec_kw={"width_ratios": [3.1, 1.25], "wspace": 0.12})
    h = 0.72
    for N in range(1, 9):
        y = 8 - N
        x = 0.0
        for l in range(0, (N - 1) // 2 + 1):
            n = N - l
            cap = 2 * (2 * l + 1)
            a.add_patch(Rectangle((x, y - h / 2), cap, h, fc=LCOL[l], ec="white", lw=1.5, alpha=0.92))
            a.text(x + cap / 2, y, f"{n}{names[l]}", ha="center", va="center", color="white",
                   fontsize=10 if cap > 2 else 9, fontweight="bold")
            x += cap
        size = ((N + 1) // 2) ** 2
        a.text(x + 0.7, y, rf"{size} harmonic{'s' if size > 1 else ''} $\times$ 2 spin $=$ {2 * size}", va="center", fontsize=9.5, color=INK)
    a.set_xlim(0, 52); a.set_ylim(-0.6, 7.6)
    a.set_yticks(range(8)); a.set_yticklabels([str(N) for N in range(8, 0, -1)])
    a.set_ylabel(r"$N=n+l$"); a.set_xlabel("states of the subshells")
    a.set_xticks([0, 2, 8, 18, 32]); a.spines["left"].set_visible(True)
    a.tick_params(axis="y", length=0)
    for l in range(4):
        a.plot([], [], "s", color=LCOL[l], label=f"$l={l}$ ({names[l]})")
    a.legend(frameon=False, loc="upper right", fontsize=9)
    periods = [2, 8, 8, 18, 18, 32, 32]
    for k, p in enumerate(periods):
        b.barh(7 - k, p, height=h, color=GREY, alpha=0.75)
        b.text(p + 0.8, 7 - k, str(p), va="center", fontsize=10)
    b.set_ylim(-0.6 + 0.0, 7.6 - 1.0 + 0.0 + 0.0)
    b.set_ylim(-0.6, 7.6)
    b.set_yticks(range(1, 8)); b.set_yticklabels([str(k) for k in range(7, 0, -1)])
    b.set_ylabel("period of the table"); b.set_xlabel("elements per period")
    b.set_xlim(0, 40); b.set_xticks([0, 10, 20, 30]); b.spines["left"].set_visible(True)
    b.tick_params(axis="y", length=0)
    a.text(0, 8.15, "(a) Madelung groups", fontsize=10.5, color=INK, fontweight="bold")
    b.text(0, 8.15, "(b) periods", fontsize=10.5, color=INK, fontweight="bold")
    save(fig, "fig3_madelung_groups")


# ---------------------------------------------------------------- 4 spectrum
def fig4():
    from p3_spindle import lowest_roots, do_eq25
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    x = np.linspace(0.6, 5.2, 300)
    cols = {1 / 3: TEAL, 1 / 2: BLUE, 2 / 3: ORANGE}
    labs = {1 / 3: r"$\mu=1/3$", 1 / 2: r"$\mu=1/2$", 2 / 3: r"$\mu=2/3$"}
    for mu, c in cols.items():
        ax.plot(x, 2 * mu**2 * (x + 1 / (2 * mu)) * (x + 1 / (2 * mu) - 1), color=c, lw=1.7, label=labs[mu])
    for mu, c in cols.items():
        for l, mk in ((0, "o"), (1, "s")):
            kappa = (2 * l + 1) / (2 * mu)
            want = [do_eq25(mu, l, nr) / mu**2 for nr in range(3)]
            got = lowest_roots(kappa, 3, want[-1] * 1.05)
            for nr, E in enumerate(got):
                Neff = nr + l + 1 + (1 / mu - 1) * l
                ax.plot(Neff, 2 * mu**2 * E, mk, color=c, mfc=c if l == 0 else "white", mec=c, ms=6.5, mew=1.5, zorder=4)
    # Madelung couplings: the mu = 1/2 curve at integer N (the x-axis ticks)
    ax.set_xticks(range(1, 6))
    ax.text(5.1, 0.5, r"$\mu=1/2$ at integer $N$:" "\n" r"$Z_NR=N(N+1)/2$", ha="right", va="bottom", fontsize=10.5, color=BLUE)
    ax.plot([], [], "o", color=INK, label=r"numerical, $l=0$")
    ax.plot([], [], "s", mfc="white", mec=INK, color=INK, mew=1.5, label=r"numerical, $l=1$")
    ax.set_xlabel(r"$N=n_r+l+1+(\mu^{-1}-1)\,l$")
    ax.set_ylabel(r"coupling $Z_N R$")
    ax.set_xlim(0.6, 5.2); ax.set_ylim(0, 17)
    ax.legend(frameon=False, loc="upper left")
    save(fig, "fig4_spectrum")


# ---------------------------------------------------------------- 5 Green's function
def fig5():
    import mpmath as mp
    import payoff_green_check as P
    R, ZR = P.R, 1.7
    Theta = np.deg2rad(70.0)
    d = np.array([np.sin(Theta), 0, np.cos(Theta)])
    x = np.array([0, 0, 1.0]) * R * 1.0
    rs = R * np.exp(np.linspace(-1.9, 1.9, 61))
    closed = np.array([float(P.G_closed(x, r * d, ZR)) for r in rs])
    far = np.where(np.abs(np.log(rs / R)) > 0.45)[0][::3]
    ref = np.array([float(P.G_radial(x, rs[i] * d, ZR)) for i in far])
    err = np.abs(closed[far] - ref) / np.abs(ref)
    print("max rel err on figure:", err.max())
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.2, 5.2), sharex=True, gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.08})
    a.plot(rs / R, np.abs(closed), color=BLUE, lw=1.9, label="closed form (Theorem D)")
    a.plot(rs[far] / R, np.abs(ref), "o", color=ORANGE, ms=5.5, mfc="white", mew=1.6, label="partial-wave sum")
    a.set_xscale("log"); a.set_yscale("log")
    from matplotlib.ticker import NullFormatter
    from matplotlib.ticker import FixedLocator, FuncFormatter
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    ylo, yhi = a.get_ylim()
    print("fig5 |G| range:", ylo, yhi)
    a.set_xticks([0.2, 0.5, 1, 2, 5]); a.xaxis.set_major_formatter(plain)
    a.yaxis.set_major_locator(FixedLocator([t for t in (0.01, 0.02, 0.03, 0.05, 0.1, 0.2) if ylo <= t <= yhi]))
    a.yaxis.set_major_formatter(plain)
    for axx in (a.xaxis, a.yaxis):
        axx.set_minor_formatter(NullFormatter())
    print("sign of G on ray:", set(np.sign(closed)))
    a.set_ylabel(r"$|G(x,y)|$"); a.legend(frameon=False)
    b.plot(rs[far] / R, np.maximum(err, 1e-17), "o", color=TEAL, ms=5)
    b.set_yscale("log"); b.set_ylim(1e-18, 1e-10); b.set_yticks([1e-17, 1e-14, 1e-11])
    b.set_ylabel("relative error"); b.set_xlabel(r"$|y|/R$ along a ray at $70^\circ$ from $x$, $|x|=R$")
    for ax in (a, b):
        ax.spines["left"].set_visible(True)
    b.axvspan(np.exp(-0.45), np.exp(0.45), color=GREY, alpha=0.15, lw=0)
    save(fig, "fig5_green")


if __name__ == "__main__":
    which = [a for a in sys.argv[1:] if a != "--check"] or ["1", "2", "3", "4", "5"]
    for w in which:
        {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5}[w]()
