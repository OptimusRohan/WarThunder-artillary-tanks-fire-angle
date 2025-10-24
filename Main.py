import math

g = 9.8  # m/s^2

def simulate_range(v0, angle_rad, mass, Cd, area, rho=1.225, dt=0.01, max_time=300.0):
    """
    Simulate projectile motion with quadratic drag using RK4.
    Returns (range_m, flight_time_s).
    """
    # Precompute drag coefficient factor C = 0.5 * rho * Cd * A
    C = 0.5 * rho * Cd * area

    # state: [x, y, vx, vy]
    vx0 = v0 * math.cos(angle_rad)
    vy0 = v0 * math.sin(angle_rad)
    s = [0.0, 0.0, vx0, vy0]

    def deriv(state):
        x, y, vx, vy = state
        v = math.hypot(vx, vy)
        # drag acceleration components (quadratic drag, opposite v)
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
        # RK4 step
        k1 = deriv(s)
        s1 = [s[i] + 0.5 * dt * k1[i] for i in range(4)]
        k2 = deriv(s1)
        s2 = [s[i] + 0.5 * dt * k2[i] for i in range(4)]
        k3 = deriv(s2)
        s3 = [s[i] + dt * k3[i] for i in range(4)]
        k4 = deriv(s3)
        s = [s[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(4)]

        t += dt
        # stop when projectile crosses ground (y < 0) after launch
        if s[1] < 0.0:
            # linear interpolation back to y=0 for better range estimate
            x_prev, y_prev = prev_state[0], prev_state[1]
            x_now, y_now = s[0], s[1]
            if y_now == y_prev:
                x_ground = x_now
            else:
                frac = -y_prev / (y_now - y_prev)
                x_ground = x_prev + frac * (x_now - x_prev)
            return x_ground, t
        prev_state = s[:]

    # if max_time reached without hitting ground, return current x
    return s[0], t

def find_angle_for_range(target_range, v0, mass, Cd, area, rho=1.225, dt=0.01, tol=1.0, max_iter=50):
    """
    Find a firing angle in degrees giving range ~ target_range (within tol meters)
    Uses bisection between 0.001 and 89.999 degrees. Returns angle_deg or None if no hit.
    """
    low_deg = 0.001
    high_deg = 89.999

    r_low, _ = simulate_range(v0, math.radians(low_deg), mass, Cd, area, rho, dt)
    r_high, _ = simulate_range(v0, math.radians(high_deg), mass, Cd, area, rho, dt)

    # If both less than target or both greater then no bracket (but we can still try)
    for i in range(max_iter):
        mid_deg = 0.5 * (low_deg + high_deg)
        r_mid, _ = simulate_range(v0, math.radians(mid_deg), mass, Cd, area, rho, dt)
        # print(i, low_deg, high_deg, r_low, r_high, mid_deg, r_mid)
        if abs(r_mid - target_range) <= tol:
            return mid_deg
        # Decide which half to keep: compare mid range to target
        # We want to bracket the root of f(angle) = range(angle) - target_range
        if (r_low - target_range) * (r_mid - target_range) <= 0:
            high_deg = mid_deg
            r_high = r_mid
        else:
            low_deg = mid_deg
            r_low = r_mid
    return None

if __name__ == "__main__":
    # Example parameters (change to your projectile's properties)
    vehicle = "m44"
    v0 = 563.0  # muzzle velocity m/s

    # Example shell geometry & mass (approximate; change to your real values)
    shell_diameter = 0.155  # meters (155 mm)
    area = math.pi * (shell_diameter / 2) ** 2
    mass = 43.1  # kg (example shell mass — please replace with actual)
    Cd = 0.295  # drag coefficient (example; shells vary)
    rho = 1.225  # air density at sea level kg/m^3

    target = float(input("Target distance (m): "))

    angle_deg = find_angle_for_range(target, v0, mass, Cd, area, rho, dt=0.01, tol=1.0)
    if angle_deg is None:
        print("No firing angle found within search limits (target may be out of range).")
    else:
        rng, flight_t = simulate_range(v0, math.radians(angle_deg), mass, Cd, area, rho, dt=0.01)
        print(f"Estimated firing angle: {angle_deg:.3f}°  -> range {rng:.1f} m, flight time {flight_t:.1f} s")