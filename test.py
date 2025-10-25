from flask import Flask, request, jsonify
from flask_cors import CORS
from datetime import datetime
import numpy as np
from astropy.time import Time
from scipy.optimize import minimize, minimize_scalar

app = Flask(__name__)
CORS(app)
class CometOrbitCalculator:'''efeed'''
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
            M_earth = 357.52911 + 35999.05029 * t - 0.0001537 * t**2
            C_earth = (1.914602 - 0.004817 * t - 0.000014 * t**2) * np.sin(np.radians(M_earth)) + \
                     (0.019993 - 0.000101 * t) * np.sin(2 * np.radians(M_earth)) + \
                     0.000289 * np.sin(3 * np.radians(M_earth))
            nu_earth = M_earth + C_earth
            R_earth = 1.000001018 * (1 - 0.01670862**2) / (1 + 0.01670862 * np.cos(np.radians(nu_earth)))
            x_earth = R_earth * np.cos(np.radians(nu_earth))
            y_earth = R_earth * np.sin(np.radians(nu_earth))
            z_earth = 0.0
            return np.array([x_earth, y_earth, z_earth])

        initial_guess = [1.52366, 0.0934, 1.85, 49.58, 286.50, Time('2025-12-01 00:00:00').jd]

        def residuals(params):
            a, e, i, Omega, omega, T = params
            total_error = 0
            for obs in self.observations:
                ra_calc, dec_calc = self.calculate_accurate_position(a, e, i, Omega, omega, T, obs['jd'], get_earth_position)
                ra_obs_deg = obs['ra'] * 15
                ra_error = min(abs(ra_calc - ra_obs_deg),
                             abs(ra_calc - ra_obs_deg + 360),
                             abs(ra_calc - ra_obs_deg - 360))
                dec_error = dec_calc - obs['dec']
                total_error += ra_error**2 + dec_error**2
            return total_error

        bounds = [
            (1.3, 1.7), (0.08, 0.11), (1.0, 3.0),
            (40.0, 60.0), (280.0, 300.0),
            (Time('2025-11-01 00:00:00').jd, Time('2025-12-31 23:59:59').jd)
        ]

        result = minimize(residuals, initial_guess, method='L-BFGS-B', bounds=bounds,
                        options={'ftol': 1e-12, 'gtol': 1e-10, 'maxiter': 1000})

        a, e, i, Omega, omega, T = result.x
        return [a, e, i, Omega, omega, T]

    def calculate_accurate_position(self, a, e, i, Omega, omega, T, jd, earth_pos_func):
        t = jd - T
        n = np.sqrt(self.GM_sun / a**3)
        M = n * t
        E = self.solve_kepler_accurate(M, e)
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2), np.sqrt(1 - e) * np.cos(E/2))
        r = a * (1 - e * np.cos(E))

        i_rad = np.radians(i)
        Omega_rad = np.radians(Omega)
        omega_rad = np.radians(omega)

        x_orb = r * np.cos(nu)
        y_orb = r * np.sin(nu)

        x_hel = (np.cos(omega_rad) * np.cos(Omega_rad) - np.sin(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * x_orb + \
                (-np.sin(omega_rad) * np.cos(Omega_rad) - np.cos(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * y_orb
        y_hel = (np.cos(omega_rad) * np.sin(Omega_rad) + np.sin(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * x_orb + \
                (-np.sin(omega_rad) * np.sin(Omega_rad) + np.cos(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * y_orb
        z_hel = (np.sin(omega_rad) * np.sin(i_rad)) * x_orb + (np.cos(omega_rad) * np.sin(i_rad)) * y_orb

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

    def solve_kepler_accurate(self, M, e, iterations=10):
        E = M
        for _ in range(iterations):
            delta_E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
            E -= delta_E
            if abs(delta_E) < 1e-12:
                break
        return E

    def calculate_earth_distance(self, a, e, i, Omega, omega, T, jd):
        def get_earth_position(jd):
            t = (jd - 2451545.0) / 36525.0
            M_earth = 357.52911 + 35999.05029 * t - 0.0001537 * t**2
            C_earth = (1.914602 - 0.004817 * t - 0.000014 * t**2) * np.sin(np.radians(M_earth)) + \
                     (0.019993 - 0.000101 * t) * np.sin(2 * np.radians(M_earth)) + \
                     0.000289 * np.sin(3 * np.radians(M_earth))
            nu_earth = M_earth + C_earth
            R_earth = 1.000001018 * (1 - 0.01670862**2) / (1 + 0.01670862 * np.cos(np.radians(nu_earth)))
            x_earth = R_earth * np.cos(np.radians(nu_earth))
            y_earth = R_earth * np.sin(np.radians(nu_earth))
            z_earth = 0.0
            return np.array([x_earth, y_earth, z_earth])

        t = jd - T
        n = np.sqrt(self.GM_sun / a**3)
        M = n * t
        E = self.solve_kepler_accurate(M, e)
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2), np.sqrt(1 - e) * np.cos(E/2))
        r = a * (1 - e * np.cos(E))

        i_rad = np.radians(i)
        Omega_rad = np.radians(Omega)
        omega_rad = np.radians(omega)

        x_orb = r * np.cos(nu)
        y_orb = r * np.sin(nu)

        x_hel = (np.cos(omega_rad) * np.cos(Omega_rad) - np.sin(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * x_orb + \
                (-np.sin(omega_rad) * np.cos(Omega_rad) - np.cos(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * y_orb
        y_hel = (np.cos(omega_rad) * np.sin(Omega_rad) + np.sin(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * x_orb + \
                (-np.sin(omega_rad) * np.sin(Omega_rad) + np.cos(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * y_orb
        z_hel = (np.sin(omega_rad) * np.sin(i_rad)) * x_orb + (np.cos(omega_rad) * np.sin(i_rad)) * y_orb

        x_earth, y_earth, z_earth = get_earth_position(jd)
        distance = np.sqrt((x_hel - x_earth)**2 + (y_hel - y_earth)**2 + (z_hel - z_earth)**2)
        return distance

def calculate_orbit_from_observations(observations_array):
    calculator = CometOrbitCalculator()
    for obs in observations_array:
        ra, dec, datetime_str = obs
        calculator.add_observation(ra, dec, datetime_str)

    orbital_elements = calculator.calculate_orbital_elements()
    a, e, i, Omega, omega, T = orbital_elements

    start_jd = Time('2025-01-01 00:00:00').jd
    end_jd = Time('2027-01-01 00:00:00').jd

    def distance_function(jd):
        return calculator.calculate_earth_distance(a, e, i, Omega, omega, T, jd)

    result = minimize_scalar(distance_function, bounds=(start_jd, end_jd), method='bounded')
    min_distance_jd = result.x
    min_distance = result.fun
    min_distance_date = Time(min_distance_jd, format='jd').datetime

    return orbital_elements, min_distance_date, min_distance

@app.route('/')
def home():
    return jsonify({
        "message": "🌠 Comet Orbit Calculation API",
        "version": "1.0",
        "endpoints": {
            "POST /api/calculate-orbit": "Calculate comet orbit from observations",
            "GET /api/example": "Get example observation data"
        }
    })

@app.route('/api/calculate-orbit', methods=['POST'])
def calculate_orbit():
    try:
        data = request.json
        if not data or 'observations' not in data:
            return jsonify({"success": False, "error": "No observations provided"}), 400

        observations = data['observations']
        if len(observations) < 3:
            return jsonify({"success": False, "error": "At least 3 observations are required"}), 400

        observations_array = []
        for obs in observations:
            if 'ra' not in obs or 'dec' not in obs or 'datetime' not in obs:
                return jsonify({
                    "success": False,
                    "error": "Each observation must contain 'ra', 'dec', and 'datetime' fields"
                }), 400
            observations_array.append([float(obs['ra']), float(obs['dec']), obs['datetime']])

        orbital_elements, min_distance_date, min_distance = calculate_orbit_from_observations(observations_array)
        a, e, i, Omega, omega, T = orbital_elements

        result = {
            "success": True,
            "orbital_elements": {
                "a": float(a), "e": float(e), "i": float(i),
                "Omega": float(Omega), "omega": float(omega), "T": float(T)
            },
            "closest_approach": {
                "date": min_distance_date.strftime("%Y-%m-%d %H:%M:%S"),
                "distance_au": float(min_distance),
                "distance_km": float(min_distance * 149597870.7)
            },
            "period": float(2 * np.pi * np.sqrt(a**3 / 0.000295912) / 365.25)
        }

        return jsonify(result)

    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500

@app.route('/api/example', methods=['GET'])
def get_example():
    example_observations = [
        [15.30977, -18.61633, "2025-10-25 00:00:00"],
        [15.40572, -18.99403, "2025-10-27 00:00:00"],
        [15.50238, -19.36158, "2025-10-29 00:00:00"],
        [15.59917, -19.71861, "2025-10-31 00:00:00"],
        [15.86444, -20.06417, "2025-11-02 00:00:00"],
    ]

    return jsonify({
        "example_observations": [
            {"ra": obs[0], "dec": obs[1], "datetime": obs[2]} for obs in example_observations
        ]
    })

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8000, debug=True)
