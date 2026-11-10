"""ViennaPS smoke test: isotropic deposition into a 2D trench, headless, no VTK window."""
import viennaps as ps
ps.setDimension(2)
ps.Logger.setLogLevel(ps.LogLevel.WARNING)
grid = 0.02   # um
domain = ps.Domain(grid, 1.0, ps.BoundaryType.REFLECTIVE_BOUNDARY)
# trench: width 0.2 um, depth 2.0 um -> aspect ratio 10
ps.MakeTrench(domain=domain, trenchWidth=0.2, trenchDepth=2.0, trenchTaperAngle=0.0, maskHeight=0.0).apply()
n0 = len(domain.getLevelSets())
model = ps.SingleParticleProcess(rate=0.01, stickingProbability=0.1, sourceExponent=1.0)
domain.duplicateTopLevelSet(ps.Material.Al2O3)
ps.Process(domain, model, 10.0).apply()   # 10 time units of growth
domain.saveSurfaceMesh("figures/viennaps_trench_smoke.vtp", True)
print("levelsets before/after:", n0, len(domain.getLevelSets()))
print("ViennaPS smoke OK")

# quick picture of the result: initial trench vs deposited surface
import numpy as np, viennals as ls
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4, 6))
for lsd, lab, c in zip(domain.getLevelSets(), ["Si trench", "Al2O3 after deposition"], ["k", "tab:orange"]):
    m = ls.Mesh(); ls.ToSurfaceMesh(lsd, m).apply()
    pts = np.array(m.getNodes())[:, :2]
    ax.plot(pts[:, 0], pts[:, 1], ".", ms=2, color=c, label=lab)
ax.set_aspect("equal"); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
ax.set_title("ViennaPS smoke: AR 10 trench, s = 0.1", fontsize=9); ax.legend(fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig("figures/viennaps_trench_smoke.png", dpi=150)
print("saved figures/viennaps_trench_smoke.png")
