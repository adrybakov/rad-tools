import radtools as rad

l = rad.lattice_example(f"BCT2")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "bct2_wigner-seitz.png",
    elev=41,
    azim=59,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=41, azim=59)
