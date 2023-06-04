import radtools as rad

l = rad.lattice_example(f"FCC")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "fcc_wigner-seitz.png",
    elev=46,
    azim=19,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=46, azim=19)
