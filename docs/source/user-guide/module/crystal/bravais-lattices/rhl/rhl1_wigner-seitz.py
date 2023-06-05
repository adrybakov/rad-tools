import radtools as rad

l = rad.lattice_example(f"RHL1")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "rhl1_wigner-seitz.png",
    elev=19,
    azim=-19,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=19, azim=-19)
