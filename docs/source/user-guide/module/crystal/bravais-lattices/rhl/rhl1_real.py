import radtools as rad

l = rad.lattice_example(f"RHL1")
l.plot("primitive")
# Save an image:
l.savefig(
    "rhl1_real.png",
    elev=35,
    azim=52,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=35, azim=52)
