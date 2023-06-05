import radtools as rad

l = rad.lattice_example(f"TRI1a")
l.plot("wigner-seitz")
# Save an image:
l.savefig(
    "tri1a_wigner-seitz.png",
    elev=9,
    azim=18,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=9, azim=18)
