import radtools as rad

l = rad.lattice_example(f"MCLC1")
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
    "mclc1_real.png",
    elev=29,
    azim=83,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=29, azim=83)
