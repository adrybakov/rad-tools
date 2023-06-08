import radtools as rad

l = rad.lattice_example(f"ORC")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "orc_wigner-seitz.png",
    elev=20,
    azim=30,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=20, azim=30)
