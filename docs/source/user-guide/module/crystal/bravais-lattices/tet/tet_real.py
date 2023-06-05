import radtools as rad

l = rad.lattice_example(f"TET")
l.plot("primitive")
# Save an image:
l.savefig(
    "tet_real.png",
    elev=30,
    azim=30,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=30, azim=30)
