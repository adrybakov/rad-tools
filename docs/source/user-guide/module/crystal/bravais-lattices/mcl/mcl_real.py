import radtools as rad

l = rad.lattice_example(f"MCL")
l.plot("primitive")
# Save an image:
l.savefig(
    "mcl_real.png",
    elev=25,
    azim=40,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=25, azim=40)
