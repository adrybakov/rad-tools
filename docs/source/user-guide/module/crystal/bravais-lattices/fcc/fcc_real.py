import radtools as rad

l = rad.lattice_example(f"FCC")
l.plot(
    "primitive",
    label="primitive",
)
l.legend()
l.plot(
    "conventional",
    label="conventional",
    colour="black"
)
l.legend()
# Save an image:
l.savefig(
    "fcc_real.png",
    elev=28,
    azim=23,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=28, azim=23)
