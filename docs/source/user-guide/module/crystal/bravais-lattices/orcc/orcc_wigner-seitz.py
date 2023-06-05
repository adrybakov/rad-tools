import radtools as rad

l = rad.lattice_example(f"ORCC")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "orcc_wigner-seitz.png",
    elev=33,
    azim=-40,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=33, azim=-40)
