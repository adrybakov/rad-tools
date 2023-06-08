import radtools as rad

l = rad.lattice_example(f"TRI2a")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "tri2a_wigner-seitz.png",
    elev=30,
    azim=62,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=30, azim=62)
