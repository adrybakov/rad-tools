import radtools as rad

l = rad.lattice_example(f"CUB")
l.plot("primitive")
# Save an image:
l.savefig(
    "cub_real.png",
    elev=28,
    azim=23,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=28, azim=23)
