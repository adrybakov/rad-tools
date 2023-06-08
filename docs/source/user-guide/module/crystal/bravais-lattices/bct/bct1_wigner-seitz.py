import radtools as rad

l = rad.lattice_example(f"BCT1")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "bct1_wigner-seitz.png",
    elev=26,
    azim=59,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=26, azim=59)
