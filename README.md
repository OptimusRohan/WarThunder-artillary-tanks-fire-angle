# WarThunder-artillery-tanks-fire-angle
This program should give the angle to fire at when using either a M44 or M55 to hit a given target.
I may later add other vehicles as well.

- I have a version that currently uses a rough estimate for air resistance and gets numbers but I have yet to test its accuracy
- Also I only have the M44 right now, I will get this one working before I start working on others
- I need to take air resistance into account to get any kind of accuracy
- I am happy to see what other people might add
- This also assumes that you are at the same height and the person you are firing at is not moving

- please report bugs here https://docs.google.com/forms/d/e/1FAIpQLScKJ0ym2Ly-CYY0b7BGp2pdREYN7-pYbjBEL_1iDBfcKzGMCA/viewform?usp=publish-editor



I worked the math and assuming no air resistance the angle of fire will be
v = fire velocity
g = gravity
x = distance from target

angle = arcsin(gx/(v^2))/2
this means half the arcsin of gravity times distance form target
divided by the fire velocity squired