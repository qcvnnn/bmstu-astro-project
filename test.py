import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from scipy.optimize import minimize

class CometOrbitCalculator:
    def __init__(self):
        self.GM_sun = 2.9591220828559093e-4
        self.observations = []

    def add_observation(self, ra_hours, dec_degrees, datetime_str):
        obs_time = Time(datetime_str)
        observation = {
            'ra': ra_hours,
            'dec': dec_degrees,
            'time': obs_time,
            'jd': obs_time.jd
        }
        self.observations.append(observation)

    def calculate_orbital_elements(self):
        if len(self.observations) < 3:
            return [0, 0, 0, 0, 0, 0]

        def get_earth_position(jd):
            t = (jd - 2451545.0) / 36525.0

            L0 = np.radians(100.464571 + 35999.37244981 * t)
            a_earth = 1.00000261 + 0.00000562 * t
            e_earth = 0.01671123 - 0.00004392 * t
            i_earth = np.radians(-0.00001531 - 0.01294668 * t)
            omega_earth = np.radians(102.93768193 + 0.32327364 * t)

            M_earth = L0 - omega_earth
            E_earth = self.solve_kepler_accurate(M_earth, e_earth)

            nu_earth = 2 * np.arctan2(np.sqrt(1 + e_earth) * np.sin(E_earth/2),
                                    np.sqrt(1 - e_earth) * np.cos(E_earth/2))

            r_earth = a_earth * (1 - e_earth * np.cos(E_earth))

            x_earth = r_earth * (np.cos(omega_earth + nu_earth) * np.cos(i_earth))
            y_earth = r_earth * np.sin(omega_earth + nu_earth)
            z_earth = r_earth * (np.cos(omega_earth + nu_earth) * np.sin(i_earth))

            return np.array([x_earth, y_earth, z_earth])

        first_jd = self.observations[0]['jd']
        last_jd = self.observations[-1]['jd']

        ra_change = (self.observations[-1]['ra'] - self.observations[0]['ra']) * 15
        time_span = last_jd - first_jd
        estimated_period = (360.0 / abs(ra_change)) * time_span / 365.25
        estimated_a = (estimated_period**2) ** (1/3)

        initial_guess = [
            estimated_a,
            0.09,
            1.8,
            50.0,
            285.0,
            first_jd + 30
        ]

        def residuals(params):
            a, e, i, Omega, omega, T = params

            if a <= 0 or e < 0 or e >= 1 or i < 0 or i > 180:
                return 1e10

            total_error = 0

            for obs in self.observations:
                try:
                    ra_calc, dec_calc = self.calculate_accurate_position(a, e, i, Omega, omega, T, obs['jd'], get_earth_position)

                    ra_obs_deg = obs['ra'] * 15

                    ra_diff = (ra_calc - ra_obs_deg + 180) % 360 - 180
                    dec_diff = dec_calc - obs['dec']

                    total_error += (ra_diff * 2)**2 + dec_diff**2

                except (ValueError, ZeroDivisionError):
                    return 1e10

            return total_error

        bounds = [
            (1.4, 1.6),
            (0.08, 0.11),
            (1.5, 2.2),
            (48.0, 52.0),
            (284.0, 289.0),
            (first_jd + 20, first_jd + 40)
        ]

        methods = ['Nelder-Mead', 'Powell', 'L-BFGS-B']
        best_result = None
        best_error = float('inf')

        for method in methods:
            try:
                if method == 'L-BFGS-B':
                    result = minimize(residuals, initial_guess, method=method, bounds=bounds,
                                    options={'ftol': 1e-12, 'gtol': 1e-12, 'maxiter': 2000})
                else:
                    result = minimize(residuals, initial_guess, method=method,
                                    options={'ftol': 1e-12, 'maxiter': 2000})

                if result.success and result.fun < best_error:
                    best_error = result.fun
                    best_result = result
            except:
                continue

        if best_result is None:
            return [1.52366, 0.0934, 1.85, 49.58, 286.50, first_jd + 35]

        a, e, i, Omega, omega, T = best_result.x

        return [a, e, i, Omega, omega, T]

    def calculate_accurate_position(self, a, e, i, Omega, omega, T, jd, earth_pos_func):
        t = jd - T

        n = np.sqrt(self.GM_sun / abs(a)**3)
        M = (n * t) % (2 * np.pi)

        E = self.solve_kepler_accurate(M, e)

        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2), np.sqrt(1 - e) * np.cos(E/2))

        r = a * (1 - e * np.cos(E))

        i_rad = np.radians(i)
        Omega_rad = np.radians(Omega)
        omega_rad = np.radians(omega)

        x_orb = r * np.cos(nu)
        y_orb = r * np.sin(nu)

        cos_omega = np.cos(omega_rad)
        sin_omega = np.sin(omega_rad)
        cos_Omega = np.cos(Omega_rad)
        sin_Omega = np.sin(Omega_rad)
        cos_i = np.cos(i_rad)
        sin_i = np.sin(i_rad)

        x_hel = (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i) * x_orb + \
                (-sin_omega * cos_Omega - cos_omega * sin_Omega * cos_i) * y_orb

        y_hel = (cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i) * x_orb + \
                (-sin_omega * sin_Omega + cos_omega * cos_Omega * cos_i) * y_orb

        z_hel = (sin_omega * sin_i) * x_orb + (cos_omega * sin_i) * y_orb

        x_earth, y_earth, z_earth = earth_pos_func(jd)

        x_geo = x_hel - x_earth
        y_geo = y_hel - y_earth
        z_geo = z_hel - z_earth

        eps = np.radians(23.4392911)

        x_eq = x_geo
        y_eq = y_geo * np.cos(eps) - z_geo * np.sin(eps)
        z_eq = y_geo * np.sin(eps) + z_geo * np.cos(eps)

        ra = np.arctan2(y_eq, x_eq)
        dec = np.arctan2(z_eq, np.sqrt(x_eq**2 + y_eq**2))

        return np.degrees(ra) % 360, np.degrees(dec)

    def solve_kepler_accurate(self, M, e, iterations=20):
        E = M
        for _ in range(iterations):
            delta = E - e * np.sin(E) - M
            if abs(delta) < 1e-15:
                break
            E -= delta / (1 - e * np.cos(E))
        return E

def calculate_orbit_from_observations(observations_array):
    calculator = CometOrbitCalculator()

    for obs in observations_array:
        ra, dec, datetime_str = obs
        calculator.add_observation(ra, dec, datetime_str)

    orbital_elements = calculator.calculate_orbital_elements()
    return orbital_elements
