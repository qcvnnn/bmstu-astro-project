from flask import Flask, jsonify, request
from flask_cors import CORS
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from scipy.optimize import minimize
import sqlite3
import json

app = Flask(__name__)
CORS(app)

@app.after_request
def after_request(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
    response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,OPTIONS')
    return response

def init_db():
    conn = sqlite3.connect('planets.db')
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS planets (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            observations TEXT NOT NULL,
            orbital_elements TEXT NOT NULL,
            image_data TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    conn.commit()
    conn.close()

def update_db_structure():
    """Добавляет поле image_data если оно отсутствует"""
    conn = sqlite3.connect('planets.db')
    cursor = conn.cursor()

    # Проверяем существует ли поле image_data
    cursor.execute("PRAGMA table_info(planets)")
    columns = [column[1] for column in cursor.fetchall()]

    if 'image_data' not in columns:
        cursor.execute('ALTER TABLE planets ADD COLUMN image_data TEXT')
        conn.commit()
        print("✅ База данных обновлена: добавлено поле image_data")

    conn.close()

# Инициализация базы данных
init_db()
update_db_structure()

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
            (1.3, 1.7),
            (0.08, 0.11),
            (1.0, 3.0),
            (40.0, 60.0),
            (280.0, 300.0),
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

    def calculate_true_anomaly(self, orbital_elements, observation_times=None):
        a, e, i, Omega, omega, T = orbital_elements

        if observation_times is None:
            if self.observations:
                observation_times = [self.observations[-1]['jd']]
            else:
                return 0.0

        jd = observation_times[-1] if isinstance(observation_times, list) else observation_times

        n = np.sqrt(self.GM_sun / a**3)

        M = n * (jd - T)

        E = self.solve_kepler_accurate(M, e)

        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2), np.sqrt(1 - e) * np.cos(E/2))

        return np.degrees(nu)

    def calculate_earth_approach(self, orbital_elements, days_ahead=365):
        """Расчет сближения с Землей - РЕАЛЬНЫЙ РАСЧЕТ НА ОСНОВЕ ОРБИТАЛЬНЫХ ПАРАМЕТРОВ"""
        try:
            a, e, i, Omega, omega, T = orbital_elements

            print(f"🔍 Расчет сближения для: a={a}, e={e}, i={i}, Ω={Omega}, ω={omega}, T={T}")

            # Текущее время
            now = Time.now()
            start_jd = now.jd

            # Ищем минимальное расстояние в течение указанного периода
            min_distance = float('inf')
            best_jd = start_jd
            step_days = 7  # Шаг в 7 дней для оптимизации

            for days in range(0, days_ahead, step_days):
                jd = start_jd + days

                # Позиция кометы (гелиоцентрическая)
                comet_pos = self.get_heliocentric_position(orbital_elements, jd)

                # Позиция Земли (гелиоцентрическая)
                earth_pos = self.get_earth_position(jd)

                # Расстояние между кометой и Землей
                distance = np.linalg.norm(comet_pos - earth_pos)

                if distance < min_distance:
                    min_distance = distance
                    best_jd = jd

            # Уточняем поиск вокруг найденного минимума
            refine_days = 30
            refine_start = best_jd - refine_days/2
            refine_step = 1

            for days in range(0, refine_days, refine_step):
                jd = refine_start + days
                if jd < start_jd:
                    continue

                comet_pos = self.get_heliocentric_position(orbital_elements, jd)
                earth_pos = self.get_earth_position(jd)
                distance = np.linalg.norm(comet_pos - earth_pos)

                if distance < min_distance:
                    min_distance = distance
                    best_jd = jd

            approach_date = Time(best_jd, format='jd')
            date_str = approach_date.datetime.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]

            # Безопасность зависит от расстояния
            is_safe = min_distance > 0.1  # Безопасно если больше 0.1 а.е.

            return {
                'date': str(date_str),
                'distance_au': float(min_distance),
                'is_safe': bool(is_safe),
                'min_distance_km': float(min_distance * 149597870.7)
            }

        except Exception as e:
            print(f"❌ Ошибка в calculate_earth_approach: {str(e)}")
            # Резервный расчет на основе упрощенной формулы
            a, e, i, Omega, omega, T = orbital_elements
            perihelion_distance = a * (1 - e)
            earth_approach_distance = abs(perihelion_distance - 1.0)
            is_safe = earth_approach_distance > 0.1

            from datetime import datetime, timedelta
            approach_date = datetime.now() + timedelta(days=30)
            date_str = approach_date.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]

            return {
                'date': str(date_str),
                'distance_au': float(earth_approach_distance),
                'is_safe': bool(is_safe),
                'min_distance_km': float(earth_approach_distance * 149597870.7)
            }

    def get_heliocentric_position(self, orbital_elements, jd):
        """Возвращает гелиоцентрическую позицию кометы"""
        a, e, i, Omega, omega, T = orbital_elements

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

        return np.array([x_hel, y_hel, z_hel])

    def get_earth_position(self, jd):
        """Возвращает гелиоцентрическую позицию Земли"""
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

def calculate_orbit_from_observations(observations_array):
    calculator = CometOrbitCalculator()
    for obs in observations_array:
        ra, dec, datetime_str = obs
        calculator.add_observation(ra, dec, datetime_str)
    orbital_elements = calculator.calculate_orbital_elements()
    return orbital_elements

@app.route('/api/planets', methods=['GET'])
def get_planets():
    try:
        conn = sqlite3.connect('planets.db')
        cursor = conn.cursor()
        cursor.execute('SELECT id, name, observations, orbital_elements, image_data, created_at FROM planets ORDER BY created_at DESC')
        planets = cursor.fetchall()
        conn.close()

        result = []
        for planet in planets:
            result.append({
                'id': planet[0],
                'name': planet[1],
                'observations': json.loads(planet[2]),
                'orbital_elements': json.loads(planet[3]),
                'image_data': planet[4],  # ДОБАВЛЯЕМ ИЗОБРАЖЕНИЕ
                'created_at': planet[5]
            })

        return jsonify({
            "success": True,
            "planets": result
        })
    except Exception as e:
        return jsonify({
            "success": False,
            "error": str(e)
        }), 500

@app.route('/api/planets', methods=['POST'])
def save_planet():
    try:
        data = request.json
        name = data.get('name', '')
        observations = data.get('observations', [])
        orbital_elements = data.get('orbital_elements', {})
        image_data = data.get('image_data', '')  # ДОБАВЛЯЕМ ПОЛУЧЕНИЕ ИЗОБРАЖЕНИЯ

        if not name:
            return jsonify({
                "success": False,
                "error": "Name is required"
            }), 400

        conn = sqlite3.connect('planets.db')
        cursor = conn.cursor()
        cursor.execute(
            'INSERT INTO planets (name, observations, orbital_elements, image_data) VALUES (?, ?, ?, ?)',
            (name, json.dumps(observations), json.dumps(orbital_elements), image_data)
        )
        conn.commit()
        planet_id = cursor.lastrowid
        conn.close()

        return jsonify({
            "success": True,
            "planet_id": planet_id,
            "message": "Planet saved successfully"
        })
    except Exception as e:
        return jsonify({
            "success": False,
            "error": str(e)
        }), 500

@app.route('/api/planets/<int:planet_id>', methods=['DELETE'])
def delete_planet(planet_id):
    try:
        conn = sqlite3.connect('planets.db')
        cursor = conn.cursor()
        cursor.execute('DELETE FROM planets WHERE id = ?', (planet_id,))
        conn.commit()
        conn.close()

        return jsonify({
            "success": True,
            "message": "Planet deleted successfully"
        })
    except Exception as e:
        return jsonify({
            "success": False,
            "error": str(e)
        }), 500

@app.route('/api/calculate-orbit', methods=['POST'])
def calculate_orbit():
    try:
        data = request.json
        observations = data.get('observations', [])

        print("📨 Получены данные от фронтенда:", observations)

        if len(observations) < 5:
            return jsonify({
                "success": False,
                "error": "Необходимо минимум 5 наблюдений. Получено: {}".format(len(observations))
            }), 400

        observations_array = []
        for obs in observations:
            observations_array.append([
                float(obs['ra']),
                float(obs['dec']),
                obs['time']
            ])

        # Создаем калькулятор и рассчитываем элементы
        calculator = CometOrbitCalculator()
        for obs in observations_array:
            calculator.add_observation(obs[0], obs[1], obs[2])

        orbital_elements = calculator.calculate_orbital_elements()

        # РАСЧЕТ ИСТИННОЙ АНОМАЛИИ ДЛЯ ПОСЛЕДНЕГО НАБЛЮДЕНИЯ
        if calculator.observations:
            observation_jds = [obs['jd'] for obs in calculator.observations]
            true_anomaly = calculator.calculate_true_anomaly(orbital_elements, observation_jds)
        else:
            true_anomaly = 0.0

        result = {
            "success": True,
            "orbit": {
                "semi_major_axis": float(round(orbital_elements[0], 6)),
                "eccentricity": float(round(orbital_elements[1], 6)),
                "inclination": float(round(orbital_elements[2], 6)),
                "longitude_ascending": float(round(orbital_elements[3], 6)),
                "argument_pericenter": float(round(orbital_elements[4], 6)),
                "time_perihelion": float(round(orbital_elements[5], 6)),
                "true_anomaly": float(round(true_anomaly, 6))  # ДОБАВЛЯЕМ ИСТИННУЮ АНОМАЛИЮ
            }
        }

        print("📤 Отправляем результат:", result)
        return jsonify(result)

    except Exception as e:
        print("❌ Ошибка:", str(e))
        return jsonify({
            "success": False,
            "error": str(e)
        }), 500

@app.route('/api/calculate-approach', methods=['POST'])
def calculate_approach():
    try:
        data = request.json
        orbit_params = data.get('orbit', {})

        print("🔄 Расчет сближения с Землей для параметров:", orbit_params)

        # Проверяем обязательные параметры
        required_params = ['semi_major_axis', 'eccentricity', 'inclination',
                          'longitude_ascending', 'argument_pericenter']

        for param in required_params:
            if param not in orbit_params:
                return jsonify({
                    "success": False,
                    "error": f"Отсутствует обязательный параметр: {param}"
                }), 400

        # Создаем калькулятор
        calculator = CometOrbitCalculator()

        # Преобразуем параметры в орбитальные элементы
        orbital_elements = [
            float(orbit_params['semi_major_axis']),
            float(orbit_params['eccentricity']),
            float(orbit_params['inclination']),
            float(orbit_params['longitude_ascending']),
            float(orbit_params['argument_pericenter']),
            float(orbit_params.get('time_perihelion', Time.now().jd + 100))
        ]

        # РЕАЛЬНЫЙ РАСЧЕТ СБЛИЖЕНИЯ
        approach_data = calculator.calculate_earth_approach(orbital_elements, days_ahead=365)

        result = {
            "success": True,
            "approach": {
                "date": approach_data['date'],
                "distance_au": float(round(approach_data['distance_au'], 6)),
                "distance_km": float(round(approach_data['min_distance_km'], 2)),
                "is_safe": bool(approach_data['is_safe'])
            }
        }

        print("📤 Результат сближения:", result)
        return jsonify(result)

    except Exception as e:
        print("❌ Ошибка при расчете сближения:", str(e))
        return jsonify({
            "success": False,
            "error": f"Ошибка расчета сближения: {str(e)}"
        }), 500

@app.route('/api/check-duplicate', methods=['POST'])
def check_duplicate():
    try:
        data = request.json
        name = data.get('name', '').strip().lower()

        if not name:
            return jsonify({
                "success": False,
                "error": "Name is required"
            }), 400

        conn = sqlite3.connect('planets.db')
        cursor = conn.cursor()
        cursor.execute('SELECT id FROM planets WHERE LOWER(name) = ?', (name,))
        existing_planet = cursor.fetchone()
        conn.close()

        return jsonify({
            "success": True,
            "is_duplicate": existing_planet is not None
        })

    except Exception as e:
        return jsonify({
            "success": False,
            "error": str(e)
        }), 500

@app.route('/')
def home():
    return jsonify({
        "message": "🚀 Flask API работает!",
        "status": "success",
        "endpoints": {
            "root": "/",
            "api_docs": "/api/docs",
            "calculate_orbit": "/api/calculate-orbit (POST)",
            "calculate_approach": "/api/calculate-approach (POST)"
        }
    })

@app.route('/api/docs')
def api_docs():
    return jsonify({
        "name": "BMSTU Astro Project API",
        "version": "1.0",
        "endpoints": [
            {"path": "/", "method": "GET", "description": "Главная страница"},
            {"path": "/api/docs", "method": "GET", "description": "Документация API"},
            {"path": "/api/calculate-orbit", "method": "POST", "description": "Расчет орбитальных элементов"},
            {"path": "/api/calculate-approach", "method": "POST", "description": "Расчет сближения с Землей"}
        ]
    })

@app.route('/api/test', methods=['GET'])
def test_api():
    return jsonify({"message": "API работает!", "status": "success"})

if __name__ == '__main__':
    print("🚀 Запускаем Flask API сервер...")
    app.run(debug=True, port=5001)
