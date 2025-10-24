# WarThunder-artillery-tanks-fire-angle
This program should give the angle to fire at when using either a M44 or M55 to hit a given target

- At the moment I am not taking air resistance into account
- I need to take air resistance into account to get any kind of accuracy
- I am happy to see what other people might add
- This also assumes that you are at the same height and the person you are firing at is not moving


I worked the math and assuming no air resistance the angle of fire will be
v = fire velocity
g = gravity
x = distance from target

angle = arcsin(gx/(v^2))/2
this means half the arcsin of gravity times distance form target
divided by the fire velocity squired