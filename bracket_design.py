"""
Parametric bracket: CadQuery CAD -> Gmsh -> CalculiX -> design optimization
=============================================================================
A wall-mounted aluminium shelf bracket is designed entirely from code:

  1. CAD: CadQuery builds the solid from 3 design parameters
     (thickness t, inner corner radius r, width w) and exports STEP/STL
  2. Mesh: Gmsh imports the STEP and builds 2nd-order tetrahedra (C3D10),
     refined at the corner where the stress peaks
  3. Solve: CalculiX static analysis (bolted to the wall, load at the tip hole)
  4. Verification: bending stress along the arm vs beam theory
  5. Stress concentration: peak stress vs fillet radius (Kt curve)
  6. Design sweep: 24 designs, lightest one that meets a safety factor of 2
     on yield and a tip deflection limit
  7. Mesh convergence of the chosen design
  8. Figures: CAD renders, 3D stress map, design space, verification

Units: mm, N, MPa, tonne. Run inside the venv with CadQuery:  python bracket_design.py
"""
import os
import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import cadquery as cq
import gmsh

# ============================ SETTINGS ============================
CCX = "ccx"
L, H = 120.0, 100.0                   # arm length, plate height [mm]
D_BOLT, D_TIP = 8.5, 8.5              # M8 clearance holes
TIP_X = L - 12.0                      # load hole position along the arm
P = 300.0                             # load at the tip hole [N] (about 30 kg)
E, NU, RHO = 69000.0, 0.33, 2.70e-9   # aluminium 6061-T6
SY = 276.0                            # yield strength [MPa]
SF_MIN = 2.0                          # required safety factor on yield
DEFL_MAX = 1.5                        # tip deflection limit [mm] (L/80, shelf bracket)
SWEEP = dict(t=[5, 6, 8], r=[3, 12, 25], w=[30, 40, 50])
N_JOBS = max(1, (os.cpu_count() or 2) - 1)
OUT = "results"
# ==================================================================

os.makedirs(OUT, exist_ok=True)
os.makedirs("runs", exist_ok=True)


# ------------------------------------------------------------------
# 1. CAD
# ------------------------------------------------------------------
def bolts(w):
    """(y, z) of the four wall bolts."""
    return [(w / 4, 30), (-w / 4, 30), (w / 4, 75), (-w / 4, 75)]


def make_bracket(t, r, w):
    """L-bracket: vertical wall plate + horizontal arm, with a smooth inner corner of radius r.
    A large radius works as a full-width gusset with no sharp re-entrant corners."""
    plate = cq.Workplane("XY").box(t, w, H, centered=(False, True, False))
    arm = cq.Workplane("XY").box(L, w, t, centered=(False, True, False)).translate((0, 0, H - t))
    body = plate.union(arm)
    if r > 0:   # concave fillet on the inner corner, where plate and arm meet
        body = body.edges(cq.selectors.NearestToPointSelector((t, 0, H - t))).fillet(r)
    for y, z in bolts(w):
        body = body.cut(cq.Workplane("YZ").workplane(offset=-1).center(y, z).circle(D_BOLT / 2).extrude(t + 2))
    body = body.cut(cq.Workplane("XY").workplane(offset=H - t - 1).center(TIP_X, 0)
                    .circle(D_TIP / 2).extrude(t + 2))
    return body


def mass_kg(body):
    return body.val().Volume() * RHO * 1e3


# ------------------------------------------------------------------
# 2. MESH (Gmsh, 2nd-order tetrahedra)
# ------------------------------------------------------------------
def mesh_step(step, t, r, w, size):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.importShapes(step)
    gmsh.model.occ.synchronize()
    f = gmsh.model.mesh.field
    # finer mesh in the corner region (fillet and gusset junction)
    box = f.add("Box")
    f.setNumber(box, "VIn", size * 0.35)
    f.setNumber(box, "VOut", size)
    f.setNumber(box, "XMin", 0)
    f.setNumber(box, "XMax", t + max(r, 3) + 6)
    f.setNumber(box, "YMin", -w)
    f.setNumber(box, "YMax", w)
    f.setNumber(box, "ZMin", H - t - max(r, 3) - 6)
    f.setNumber(box, "ZMax", H)
    f.setNumber(box, "Thickness", 6)
    f.setAsBackgroundMesh(box)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = dict(zip(tags.astype(int), coords.reshape(-1, 3)))
    _, etags, enodes = gmsh.model.mesh.getElements(3, -1)
    gmsh.finalize()
    # Gmsh 10-node tet -> CalculiX C3D10: last two mid-side nodes swap places
    conn = np.array(enodes[0], dtype=int).reshape(-1, 10)[:, [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]]
    return nodes, conn


# ------------------------------------------------------------------
# 3. CALCULIX
# ------------------------------------------------------------------
def hole_nodes(nodes, axis, center, radius, span):
    """Nodes on a cylindrical hole surface."""
    ids = []
    for n, c in nodes.items():
        if axis == "x":
            d, s = np.hypot(c[1] - center[0], c[2] - center[1]), c[0]
        else:
            d, s = np.hypot(c[0] - center[0], c[1] - center[1]), c[2]
        if abs(d - radius) < 0.02 and span[0] - 1e-6 <= s <= span[1] + 1e-6:
            ids.append(n)
    return ids


def write_inp(path, nodes, conn, t, w):
    fix = []
    for y, z in bolts(w):
        fix += hole_nodes(nodes, "x", (y, z), D_BOLT / 2, (0, t))
    tip = hole_nodes(nodes, "z", (TIP_X, 0), D_TIP / 2, (H - t, H))
    with open(path, "w") as f:
        f.write("*NODE, NSET=NALL\n")
        for n, c in nodes.items():
            f.write(f"{n}, {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}\n")
        f.write("*ELEMENT, TYPE=C3D10, ELSET=EALL\n")
        for i, e in enumerate(conn, 1):
            f.write(f"{i}, " + ", ".join(map(str, e[:8])) + ",\n" + ", ".join(map(str, e[8:])) + "\n")
        for name, ids in (("FIX", fix), ("TIP", tip)):
            f.write(f"*NSET, NSET={name}\n")
            for i in range(0, len(ids), 12):
                f.write(", ".join(map(str, ids[i:i + 12])) + ",\n")
        f.write(f"*MATERIAL, NAME=AL6061\n*ELASTIC\n{E}, {NU}\n*DENSITY\n{RHO}\n")
        f.write("*SOLID SECTION, ELSET=EALL, MATERIAL=AL6061\n")
        f.write("*BOUNDARY\nFIX, 1, 3\n*STEP\n*STATIC\n*CLOAD\n")
        for n in tip:
            f.write(f"{n}, 3, {-P / len(tip):.8f}\n")
        f.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")
    return tip


def run_ccx(job):
    subprocess.run([CCX, job], capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="1"))
    return job


def read_frd(path):
    """Nodal displacements and stress tensors from the .frd file (fixed-width format)."""
    disp, stress, block = {}, {}, None
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                block = line.split()[1]
                continue
            if line.startswith(" -3"):
                block = None
                continue
            if block in ("DISP", "STRESS") and line.startswith(" -1"):
                n = int(line[3:13])
                vals = [float(line[13 + 12 * k:25 + 12 * k]) for k in range(6 if block == "STRESS" else 3)]
                (disp if block == "DISP" else stress)[n] = np.array(vals)
    return disp, stress


def von_mises(s):
    sx, sy, sz, txy, tyz, tzx = s
    return np.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (txy ** 2 + tyz ** 2 + tzx ** 2))


def away_from_bc(c, w):
    """True if a node is far enough from bolts and load hole (Saint-Venant)."""
    for y, z in bolts(w):
        if c[0] < 20 and np.hypot(c[1] - y, c[2] - z) < 1.5 * D_BOLT:
            return False
    return np.hypot(c[0] - TIP_X, c[1]) > 1.5 * D_TIP


# ------------------------------------------------------------------
# DESIGN EVALUATION
# ------------------------------------------------------------------
def prepare(t, r, w, size=None, tag=""):
    """CAD + mesh + input deck (serial: CadQuery and Gmsh are not thread-safe)."""
    size = size or 0.9 * t
    name = f"runs/t{t}_r{r}_w{w}{tag}"
    body = make_bracket(t, r, w)
    cq.exporters.export(body, name + ".step")
    nodes, conn = mesh_step(name + ".step", t, r, w, size)
    tip = write_inp(name + ".inp", nodes, conn, t, w)
    return dict(name=name, t=t, r=r, w=w, size=size, body=body, nodes=nodes, conn=conn, tip=tip,
                mass=mass_kg(body))


def evaluate(d):
    disp, stress = read_frd(d["name"] + ".frd")
    if not stress:
        raise RuntimeError(f"CalculiX produced no results for {d['name']} (see {d['name']}.sta)")
    vm = {n: von_mises(s) for n, s in stress.items()}
    ok = {n: v for n, v in vm.items() if away_from_bc(d["nodes"][n], d["w"])}
    n_max = max(ok, key=ok.get)
    d.update(disp=disp, stress=stress, vm=vm, s_max=ok[n_max], at=d["nodes"][n_max],
             defl=-np.mean([disp[n][2] for n in d["tip"]]), elements=len(d["conn"]))
    d["SF"] = SY / d["s_max"]
    d["feasible"] = d["SF"] >= SF_MIN and d["defl"] <= DEFL_MAX
    return d


def solve_all(designs):
    with ThreadPoolExecutor(max_workers=N_JOBS) as pool:
        list(pool.map(run_ccx, [d["name"] for d in designs]))
    return [evaluate(d) for d in designs]


# ------------------------------------------------------------------
# PLOTTING HELPERS
# ------------------------------------------------------------------
def draw_cad(ax, body, color="#9fb4c7"):
    verts, tris = body.val().tessellate(0.2, 0.3)
    v = np.array([(p.x, p.y, p.z) for p in verts])
    polys = [v[list(tri)] for tri in tris]
    normals = [np.cross(p[1] - p[0], p[2] - p[0]) for p in polys]
    light = np.array([0.4, -0.6, 0.7]) / np.linalg.norm([0.4, -0.6, 0.7])
    shade = [0.45 + 0.55 * abs(np.dot(n / (np.linalg.norm(n) + 1e-12), light)) for n in normals]
    base = np.array(colors.to_rgb(color))
    ax.add_collection3d(Poly3DCollection(polys, facecolors=[base * s for s in shade], edgecolors="none"))
    ax.set_xlim(-10, L + 5)
    ax.set_ylim(-60, 60)
    ax.set_zlim(-5, H + 10)
    ax.set_box_aspect((L + 15, 120, H + 15))
    ax.view_init(elev=18, azim=-58)
    ax.set_axis_off()


TET_FACES = [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)]


def surface_tris(conn):
    count = {}
    for e in conn:
        for f in TET_FACES:
            tri = tuple(e[i] for i in f)
            key = tuple(sorted(tri))
            count[key] = count.get(key, (0, tri))
            count[key] = (count[key][0] + 1, tri)
    return [v[1] for v in count.values() if v[0] == 1]


def draw_stress(ax, d, scale, vmax):
    tris = surface_tris(d["conn"])
    norm = colors.Normalize(0, vmax)
    polys, cols = [], []
    for tri in tris:
        polys.append([d["nodes"][n] + scale * d["disp"][n] for n in tri])
        cols.append(cm.turbo(norm(np.mean([d["vm"][n] for n in tri]))))
    ax.add_collection3d(Poly3DCollection(polys, facecolors=cols, edgecolors="none"))
    ax.set_xlim(-10, L + 5)
    ax.set_ylim(-60, 60)
    ax.set_zlim(-5, H + 10)
    ax.set_box_aspect((L + 15, 120, H + 15))
    ax.view_init(elev=18, azim=-58)
    ax.set_axis_off()
    return norm


# ==================================================================
def main():
    report = []

    def log(text=""):
        print(text)
        report.append(text)

    log("PARAMETRIC ALUMINIUM BRACKET: CadQuery + Gmsh + CalculiX")
    log(f"Arm {L:.0f} mm, plate {H:.0f} mm | load {P:.0f} N at the tip hole | "
        f"6061-T6 (Sy = {SY:.0f} MPa) | targets: SF >= {SF_MIN}, tip deflection <= {DEFL_MAX} mm\n")

    # ---------------- Design sweep ----------------
    combos = list(itertools.product(SWEEP["t"], SWEEP["r"], SWEEP["w"]))
    print(f"Building {len(combos)} designs (CAD + mesh)...")
    designs = [prepare(t, r, w) for t, r, w in combos]
    print(f"Solving {len(designs)} CalculiX models, {N_JOBS} at a time...")
    designs = solve_all(designs)

    log(f"{'t':>3} {'r':>3} {'w':>3} {'mass [g]':>9} {'elements':>9} {'sigma max':>10} {'SF':>6} {'defl [mm]':>10}  ok")
    for d in designs:
        log(f"{d['t']:>3} {d['r']:>3} {d['w']:>3} {d['mass'] * 1e3:>9.1f} {d['elements']:>9} {d['s_max']:>10.1f} "
            f"{d['SF']:>6.2f} {d['defl']:>10.3f}  {'yes' if d['feasible'] else '-'}")
    feas = [d for d in designs if d["feasible"]]
    best = min(feas, key=lambda d: d["mass"])
    heavy = max(designs, key=lambda d: d["mass"])
    log(f"\nLightest design that meets both targets: t = {best['t']} mm, corner radius {best['r']} mm, width {best['w']} mm")
    log(f"  mass {best['mass'] * 1e3:.1f} g, peak stress {best['s_max']:.1f} MPa (SF {best['SF']:.2f}), "
        f"tip deflection {best['defl']:.3f} mm")
    log(f"  peak stress located at x = {best['at'][0]:.1f}, y = {best['at'][1]:.1f}, z = {best['at'][2]:.1f} mm")
    small_r = [d for d in feas if d["r"] == min(SWEEP["r"])]
    if small_r:
        lg = min(small_r, key=lambda d: d["mass"])
        log(f"  best design with the smallest corner radius ({lg['r']} mm): t = {lg['t']} mm, w = {lg['w']} mm, "
            f"{lg['mass'] * 1e3:.1f} g -> the large radius saves {(1 - best['mass'] / lg['mass']) * 100:.0f}% of the mass")

    # ---------------- Verification: bending stress along the arm ----------------
    print("\nVerification and fillet study...")
    t_v = 6
    w_v = 40
    ver = prepare(t_v, 2, w_v, tag="_verif")
    fillets = [prepare(t_v, r, w_v, tag="_kt") for r in (0.5, 1, 2, 3, 5, 8)]
    conv = [prepare(best["t"], best["r"], best["w"], size=s, tag=f"_m{k}")
            for k, s in enumerate([1.6 * best["t"], 1.1 * best["t"], 0.8 * best["t"], 0.6 * best["t"]])]
    ver, *rest = solve_all([ver] + fillets + conv)
    fillets, conv = rest[:6], rest[6:]

    I = w_v * t_v ** 3 / 12
    xs, sim = [], []
    for n, c in ver["nodes"].items():
        if abs(c[2] - H) < 1e-6 and abs(c[1]) < 0.25 * w_v and 3 * t_v < c[0] < TIP_X - 12:
            xs.append(c[0])
            sim.append(ver["stress"][n][0])
    xs, sim = np.array(xs), np.array(sim)
    theory = lambda x: P * (TIP_X - x) * (t_v / 2) / I
    err = np.mean(np.abs(sim - theory(xs)) / theory(xs)) * 100
    log(f"\nVERIFICATION (t = {t_v} mm, w = {w_v} mm, r = 2 mm): bending stress on the top of the arm vs beam theory "
        f"sigma = M c / I -> mean difference {err:.1f}%")

    s_nom = P * (TIP_X - t_v) * (t_v / 2) / I
    log(f"\nSTRESS CONCENTRATION at the inner corner (t = {t_v} mm, nominal root stress {s_nom:.1f} MPa)")
    for d in fillets:
        log(f"  fillet r = {d['r']:>4} mm: peak {d['s_max']:6.1f} MPa -> Kt = {d['s_max'] / s_nom:.2f}")

    log(f"\nMESH CONVERGENCE (chosen design)")
    for d in conv:
        log(f"  size {d['size']:5.2f} mm, {d['elements']:>7} elements: peak {d['s_max']:6.1f} MPa, "
            f"deflection {d['defl']:.4f} mm")
    ch = (conv[-1]["s_max"] - conv[-2]["s_max"]) / conv[-1]["s_max"] * 100
    log(f"  last refinement changes the peak stress by {ch:+.2f}%")

    # ---------------- Exports ----------------
    cq.exporters.export(best["body"], f"{OUT}/best_bracket.step")
    cq.exporters.export(best["body"], f"{OUT}/best_bracket.stl")

    # ================= FIGURES =================
    # 1. CAD gallery
    picks = [heavy, min(designs, key=lambda d: d["mass"]), best]
    titles = ["Heaviest design", "Lightest design", "Optimum"]
    fig = plt.figure(figsize=(15, 5))
    for k, (d, ttl) in enumerate(zip(picks, titles)):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        draw_cad(ax, d["body"], color="#e0a96d" if d is best else "#9fb4c7")
        ax.set_title(f"{ttl}: t = {d['t']}, r = {d['r']}, w = {d['w']} mm\n{d['mass'] * 1e3:.0f} g, "
                     f"SF {d['SF']:.2f}, deflection {d['defl']:.2f} mm", fontsize=10)
    fig.suptitle("Bracket geometries generated by CadQuery", fontsize=13)
    fig.tight_layout()
    fig.savefig(f"{OUT}/cad_designs.png", dpi=140)
    plt.close(fig)

    # 2. Stress map of the optimum
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(projection="3d")
    norm = draw_stress(ax, best, scale=30, vmax=best["s_max"])
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.turbo), ax=ax, shrink=0.6,
                 label="von Mises stress [MPa]")
    ax.set_title(f"Optimum design under {P:.0f} N (deformation x30)\npeak {best['s_max']:.0f} MPa, "
                 f"SF {best['SF']:.2f}, {best['elements']} C3D10 elements", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{OUT}/stress_optimum.png", dpi=150)
    plt.close(fig)

    # 3. Design space
    fig, ax = plt.subplots(figsize=(9, 5.5))
    mk = dict(zip(SWEEP["w"], ["o", "s", "^"]))
    for d in designs:
        ax.scatter(d["mass"] * 1e3, d["SF"], marker=mk[d["w"]], s=30 + 6 * d["r"],
                   c=[cm.viridis(min(d["defl"] / 2.0, 1))], edgecolors="k" if d["feasible"] else "none")
    ax.scatter(best["mass"] * 1e3, best["SF"], s=350, facecolors="none", edgecolors="tab:red", linewidths=2)
    ax.annotate("optimum", (best["mass"] * 1e3, best["SF"]), textcoords="offset points", xytext=(12, -14),
                color="tab:red")
    ax.axhline(SF_MIN, color="tab:red", ls="--", lw=1)
    ax.text(ax.get_xlim()[1], SF_MIN, f" SF = {SF_MIN}", va="center", color="tab:red", fontsize=8)
    for w, m in mk.items():
        ax.scatter([], [], marker=m, color="gray", label=f"width {w} mm")
    fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0, 2.0), cmap=cm.viridis), ax=ax,
                 label="tip deflection [mm]")
    ax.set_xlabel("Mass [g]")
    ax.set_ylabel("Safety factor on yield")
    ax.set_title(f"Design space: {len(designs)} brackets (black edge = meets SF and deflection targets; "
                 "marker size = corner radius)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(f"{OUT}/design_space.png", dpi=150)
    plt.close(fig)

    # 4. Verification and Kt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    o = np.argsort(xs)
    ax.plot(xs[o], sim[o], ".", ms=4, color="tab:red", label="FEA, top surface of the arm")
    xx = np.linspace(xs.min(), xs.max(), 50)
    ax.plot(xx, theory(xx), "k-", label="Beam theory  σ = M·c / I")
    ax.set_xlabel("Position along the arm x [mm]")
    ax.set_ylabel("Bending stress σxx [MPa]")
    ax.set_title(f"Verification (t = {t_v} mm, w = {w_v} mm): mean difference {err:.1f}%")
    ax.grid(alpha=0.3)
    ax.legend()
    ax = axes[1]
    rs = [d["r"] for d in fillets]
    kts = [d["s_max"] / s_nom for d in fillets]
    ax.plot(rs, kts, "o-", color="k")
    ax.set_xlabel("Inner fillet radius r [mm]")
    ax.set_ylabel("Stress concentration factor Kt = σ_peak / σ_nominal")
    ax.set_title("A bigger fillet lowers the peak stress at the corner")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/verification_kt.png", dpi=150)
    plt.close(fig)

    # 5. Mesh convergence
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot([d["elements"] for d in conv], [d["s_max"] for d in conv], "o-", color="k")
    ax.set_xscale("log")
    ax.set_xlabel("Number of C3D10 elements")
    ax.set_ylabel("Peak von Mises stress [MPa]")
    ax.set_title("Mesh convergence of the peak stress (optimum design)")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(f"{OUT}/mesh_convergence.png", dpi=150)
    plt.close(fig)

    with open(f"{OUT}/report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report))
    print(f"\nFigures, STEP/STL of the optimum and the report saved in '{OUT}'.")


if __name__ == "__main__":
    main()
