import radtools as rad

l = rad.lattice_example(f"BCC")
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
    "bcc_real.png",
    elev=30,
    azim=-180,
    dpi=300,
   bbox_inches="tight",
)
# Interactive plot:
l.show(elev=30, azim=-180)
