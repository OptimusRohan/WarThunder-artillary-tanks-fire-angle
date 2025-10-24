import math

vehicle = input("Are you using a (1)M44 or (2)M55? ")
while 1:
    vehicle = "M44"  # example; replace with your variable/source
    fireDistance = float(input("How far is the target (meters): "))
    gravity = 9.8

    if vehicle.lower() == "m44":
        MuzzleVelocity = 563.0  # m/s
        ratio = gravity * fireDistance / (MuzzleVelocity ** 2)

        # check domain
        if ratio < -1 or ratio > 1:
            print("No solution: required sin(2θ) = {:.6f} is outside [-1, 1].".format(ratio))
        else:
            # 2θ = asin(ratio)  -> θ = 0.5 * asin(ratio)
            theta_rad = 0.5 * math.asin(ratio)
            theta_deg = math.degrees(theta_rad)

            # complementary solution (the other angle that gives same range)
            alt_theta_deg = 90.0 - theta_deg

            print("sin(2θ) = {:.6f}".format(ratio))
            print("One firing angle = {:.6f} rad = {:.6f}°".format(theta_rad, theta_deg))
            print("Other possible firing angle = {:.6f}° (use whichever is appropriate)".format(alt_theta_deg))
