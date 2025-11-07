import math

g = 9.8  # m/s^2

def get_cd_for_speed(v, speed_of_sound, cd_subsonic, cd_transonic, cd_supersonic,
                     mach_lower=0.95, mach_upper=1.05):
    if v <= 0.0:
        return cd_subsonic
    mach = v / speed_of_sound
    if mach < mach_lower:
        return cd_subsonic
    if mach_lower <= mach <= 1.0:
        t = (mach - mach_lower) / (1.0 - mach_lower) if (1.0 - mach_lower) > 0 else 0.0
        return cd_subsonic + t * (cd_transonic - cd_subsonic)
    if 1.0 < mach <= mach_upper:
        t = (mach - 1.0) / (mach_upper - 1.0) if (mach_upper - 1.0) > 0 else 0.0
        return cd_transonic + t * (cd_supersonic - cd_transonic)
    return cd_supersonic

def simulate_range(v0, angle_rad, mass, cd_tuple, area, rho=1.225,
                   speed_of_sound=343.0, dt=0.01, max_time=300.0,
                   mach_lower=0.95, mach_upper=1.05):
    cd_sub, cd_trans, cd_sup = cd_tuple
    vx0 = v0 * math.cos(angle_rad)
    vy0 = v0 * math.sin(angle_rad)
    s = [0.0, 0.0, vx0, vy0]

    def deriv(state):
        x, y, vx, vy = state
        v = math.hypot(vx, vy)
        Cd = get_cd_for_speed(v, speed_of_sound, cd_sub, cd_trans, cd_sup, mach_lower, mach_upper)
        C = 0.5 * rho * Cd * area
        if v == 0.0:
            ax = 0.0
            ay = -g
        else:
            ax = - (C / mass) * v * vx
            ay = -g - (C / mass) * v * vy
        return [vx, vy, ax, ay]

    t = 0.0
    prev_state = s[:]
    while t < max_time:
        k1 = deriv(s)
        s1 = [s[i] + 0.5 * dt * k1[i] for i in range(4)]
        k2 = deriv(s1)
        s2 = [s[i] + 0.5 * dt * k2[i] for i in range(4)]
        k3 = deriv(s2)
        s3 = [s[i] + dt * k3[i] for i in range(4)]
        k4 = deriv(s3)
        s = [s[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(4)]

        t += dt
        if s[1] < 0.0:
            x_prev, y_prev = prev_state[0], prev_state[1]
            x_now, y_now = s[0], s[1]
            if y_now == y_prev:
                x_ground = x_now
            else:
                frac = -y_prev / (y_now - y_prev)
                x_ground = x_prev + frac * (x_now - x_prev)
            return x_ground, t
        prev_state = s[:]
    return s[0], t

def _bisect_angle_for_target(target_range, v0, mass, cd_tuple, area, rho,
                             speed_of_sound, dt, mach_lower, mach_upper,
                             low_deg, high_deg, tol, max_iter=50):
    r_low, _ = simulate_range(v0, math.radians(low_deg), mass, cd_tuple, area, rho,
                              speed_of_sound, dt, mach_lower, mach_upper)
    r_high, _ = simulate_range(v0, math.radians(high_deg), mass, cd_tuple, area, rho,
                               speed_of_sound, dt, mach_lower, mach_upper)
    # If both are extremely close already, return the midpoint
    if abs(r_low - target_range) <= tol:
        return low_deg
    if abs(r_high - target_range) <= tol:
        return high_deg

    for _ in range(max_iter):
        mid_deg = 0.5 * (low_deg + high_deg)
        r_mid, _ = simulate_range(v0, math.radians(mid_deg), mass, cd_tuple, area, rho,
                                  speed_of_sound, dt, mach_lower, mach_upper)
        if abs(r_mid - target_range) <= tol:
            return mid_deg
        # Keep the subinterval that contains a sign change for f(angle) = range(angle)-target
        if (r_low - target_range) * (r_mid - target_range) <= 0:
            high_deg = mid_deg
            r_high = r_mid
        else:
            low_deg = mid_deg
            r_low = r_mid
    return 0.5 * (low_deg + high_deg)

def find_angles_for_range(target_range, v0, mass, cd_tuple, area, rho=1.225,
                          speed_of_sound=343.0, dt=0.01, tol=1.0,
                          mach_lower=0.95, mach_upper=1.05,
                          angle_step=0.5, max_bisect_iter=50):
    """
    Find all firing angles (degrees) that give approx the target_range.
    Scans angles from 0.001 to 89.999 degrees with step=angle_step (deg),
    detects sign changes of f(angle) = simulate_range(angle) - target_range,
    and refines each bracket with bisection to ~tol meters.

    Returns a sorted list of angle (deg). If none found, returns [].
    """
    angles_found = []
    a_min = 0.001
    a_max = 89.999

    # helper to evaluate f(angle)
    def f_deg(angle_deg):
        r, _ = simulate_range(v0, math.radians(angle_deg), mass, cd_tuple, area, rho,
                              speed_of_sound, dt, mach_lower, mach_upper)
        return r - target_range

    # sample to find brackets
    sample_angles = []
    a = a_min
    while a <= a_max + 1e-9:
        sample_angles.append(a)
        a += angle_step
    if sample_angles[-1] < a_max:
        sample_angles.append(a_max)

    prev_a = sample_angles[0]
    prev_f = f_deg(prev_a)
    # check if exact at first sample
    if abs(prev_f) <= tol:
        angles_found.append(prev_a)

    for a in sample_angles[1:]:
        cur_f = f_deg(a)
        if abs(cur_f) <= tol:
            angles_found.append(a)
        elif prev_f * cur_f < 0:
            # sign change bracket found -> refine
            refined = _bisect_angle_for_target(target_range, v0, mass, cd_tuple, area, rho,
                                               speed_of_sound, dt, mach_lower, mach_upper,
                                               prev_a, a, tol, max_iter=max_bisect_iter)
            angles_found.append(refined)
        prev_a, prev_f = a, cur_f

    # remove duplicates (angles very close) and sort
    angles_unique = []
    angles_found.sort()
    for ang in angles_found:
        if not any(abs(ang - u) < 1e-3 for u in angles_unique):
            angles_unique.append(ang)
    return angles_unique

def find_high_angle_for_range(target_range, v0, mass, cd_tuple, area, rho=1.225,
                              speed_of_sound=343.0, dt=0.01, tol=1.0,
                              mach_lower=0.95, mach_upper=1.05,
                              angle_step=0.5, max_bisect_iter=50):
    """
    Convenience wrapper: returns the larger firing angle (deg) that hits target_range,
    or None if there is no solution.
    """
    angles = find_angles_for_range(target_range, v0, mass, cd_tuple, area, rho,
                                   speed_of_sound, dt, tol, mach_lower, mach_upper,
                                   angle_step, max_bisect_iter)
    if not angles:
        return None
    return max(angles)

if __name__ == "__main__":
    # Example parameters (change to your projectile's properties)
    vehicle = "m44"
    M44_v0 = 563.0  # muzzle velocity m/s for the m44
    M44_shell_diameter = 0.155  # meters (155 mm)
    area = math.pi * (M44_shell_diameter / 2) ** 2
    M44_mass = 43.1  # kg M44 shell mass
    rho = 1.225  # air density at sea level kg/m^3

    M44_cd_subsonic = 0.13
    M44_cd_transonic = 0.4
    M44_cd_supersonic = 0.25
    M44_cd_tuple = (M44_cd_subsonic, M44_cd_transonic, M44_cd_supersonic)
    speed_of_sound = 343.0  # m/s at ~20°C sea level

    target = float(input("Target distance (m): "))

    # find the high-angle solution
    high_angle = find_high_angle_for_range(target, M44_v0, M44_mass, M44_cd_tuple, area, rho,
                                           speed_of_sound=speed_of_sound, dt=0.01, tol=1.0,
                                           angle_step=0.5)
    if high_angle is None:
        print("No firing angle found within search limits (target may be out of range).")
    else:
        rng, flight_t = simulate_range(M44_v0, math.radians(high_angle), M44_mass, M44_cd_tuple, area,
                                       rho, speed_of_sound=speed_of_sound, dt=0.01)
        print(f"High firing angle: {high_angle:.3f}°  -> range {rng:.1f} m, flight time {flight_t:.1f} s")