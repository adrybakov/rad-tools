import radtools as rad

l = rad.lattice_example(f"MCLC2")
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
    "mclc2_real.png",
    elev=30,
    azim=62,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=30, azim=62)
