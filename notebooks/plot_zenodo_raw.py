"""Week 1: plot the Zenodo 3974438 PillarHall profiles as received, depth re-zeroed at the cavity entrance
(the authors' own shifted X column). Rejected fits (thickness > 100 nm or RMSE > 2) are drawn hollow."""
import sys
sys.path = [p for p in sys.path if p not in ("", ".")]
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from aldconf.data import load_pillarhall, default_path, GAP_NM

df = load_pillarhall(default_path())
print(df[df["good"]].groupby("sheet").agg(n=("thickness_nm", "size"), t_max=("thickness_nm", "max"), t_min=("thickness_nm", "min"),
                              depth_min=("depth_um", "min"), depth_max=("depth_um", "max")).to_string())
out = Path(__file__).resolve().parents[1] / "figures"
out.mkdir(exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for sh, g in df.groupby("sheet"):
    # points past the channel end sit on the open membrane-free field, not in the cavity
    inside = g["depth_um"] <= g["length_mm"] * 1000.0
    g = g[inside]
    good = g["good"]
    axes[0].plot(g.loc[good, "depth_um"], g.loc[good, "thickness_nm"], ".", ms=3, label=sh)
    axes[0].plot(g.loc[~good, "depth_um"], g.loc[~good, "thickness_nm"], "o", mfc="none", ms=4, color="k", alpha=0.4)
    ar = g["depth_um"] * 1000.0 / GAP_NM
    axes[1].plot(ar[good], g.loc[good, "thickness_nm"], ".", ms=3)
axes[0].set_xlabel("depth from cavity entrance (um)")
axes[0].set_ylabel("Al2O3 thickness (nm)")
axes[0].set_ylim(-2, 70)
axes[0].set_title("Zenodo 3974438: TMA/H2O 300 C, 500 cycles, 0.1 s pulses, 500 nm gap")
axes[0].legend(fontsize=6, ncol=2)
axes[1].set_xlabel("depth / gap (aspect ratio reached)")
axes[1].set_ylabel("Al2O3 thickness (nm)")
axes[1].set_xlim(left=-5)
axes[1].set_ylim(-2, 60)
fig.tight_layout()
fig.savefig(out / "zenodo_raw.png", dpi=150)
print("saved", out / "zenodo_raw.png")
