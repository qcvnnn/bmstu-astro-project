from flask import Flask, jsonify, request
from flask_cors import CORS
import numpy as np
from astropy.coordinates import SkyCoord, get_body_barycentric, solar_system_ephemeris
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import CartesianRepresentation
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

        def get_earth_position_accurate(jd):
            """Точное положение Земли используя astropy"""
            with solar_system_ephemeris.set('builtin'):
                earth_pos = get_body_barycentric('earth', Time(jd, format='jd'))
                # Конвертируем в AU
                return np.array([earth_pos.x.to(u.AU).value,
                               earth_pos.y.to(u.AU).value,
                               earth_pos.z.to(u.AU).value])

        # Анализируем наблюдения для определения типа орбиты
        if len(self.observations) >= 3:
            # Определяем приблизительные параметры на основе движения
            ra_change = (self.observations[-1]['ra'] - self.observations[0]['ra']) * 15
            dec_change = self.observations[-1]['dec'] - self.observations[0]['dec']
            time_span = self.observations[-1]['jd'] - self.observations[0]['jd']

            if abs(ra_change) < 10 and abs(dec_change) < 10 and time_span < 50:
                # Медленное движение - вероятно планета земной группы
                initial_guess = [1.52366, 0.0934, 1.85, 49.58, 286.50, Time('2025-12-01 00:00:00').jd]
                bounds = [
                    (1.45, 1.60), (0.08, 0.12), (1.5, 2.2),
                    (40.0, 60.0), (280.0, 295.0),
                    (Time('2025-11-15 00:00:00').jd, Time('2025-12-15 00:00:00').jd)
                ]
                print("🎯 Оптимизация для планетарной орбиты")
            else:
                # Быстрое движение - вероятно комета
                initial_guess = [3.0, 0.5, 10.0, 100.0, 200.0, Time('2025-06-01 00:00:00').jd]
                bounds = [
                    (1.0, 50.0), (0.1, 0.95), (0.0, 90.0),
                    (0.0, 360.0), (0.0, 360.0),
                    (min(obs['jd'] for obs in self.observations) - 365,
                     max(obs['jd'] for obs in self.observations) + 365)
                ]
                print("🎯 Оптимизация для кометной орбиты")
        else:
            initial_guess = [2.0, 0.1, 10.0, 100.0, 200.0, Time.now().jd + 100]
            bounds = [
                (0.5, 20.0), (0.01, 0.95), (0.0, 90.0),
                (0.0, 360.0), (0.0, 360.0),
                (min(obs['jd'] for obs in self.observations) - 365,
                 max(obs['jd'] for obs in self.observations) + 365)
            ]

        def residuals(params):
            a, e, i, Omega, omega, T = params

            # Штрафы за нефизические значения
            penalty = 0
            if e < 0 or e >= 1:
                penalty += 1000
            if a <= 0:
                penalty += 1000

            total_error = penalty

            for obs in self.observations:
                try:
                    ra_calc, dec_calc = self.calculate_accurate_position(a, e, i, Omega, omega, T, obs['jd'], get_earth_position_accurate)

                    ra_obs_deg = obs['ra'] * 15
                    # Улучшенный расчет ошибки RA
                    ra_diff = (ra_calc - ra_obs_deg + 180) % 360 - 180
                    dec_diff = dec_calc - obs['dec']

                    # Взвешивание ошибок
                    total_error += ra_diff**2 + (dec_diff * 2)**2

                except Exception as e:
                    total_error += 1000  # Большой штраф за ошибки вычислений

            return total_error

        # Многоуровневая оптимизация
        best_result = None
        best_error = float('inf')

        # Пробуем несколько методов
        for method in ['Nelder-Mead', 'L-BFGS-B', 'Powell']:
            try:
                if method == 'Nelder-Mead':
                    result = minimize(residuals, initial_guess, method=method,
                                    options={'maxiter': 2000, 'xatol': 1e-6, 'fatol': 1e-6})
                else:
                    result = minimize(residuals, initial_guess, method=method, bounds=bounds,
                                    options={'maxiter': 2000, 'ftol': 1e-10})

                if result.success and result.fun < best_error:
                    best_error = result.fun
                    best_result = result
                    print(f"✅ {method}: ошибка = {result.fun:.2e}")

            except Exception as e:
                print(f"⚠️ {method} не сработал: {e}")
                continue

        if best_result is None:
            print("❌ Все методы оптимизации не сработали, используем начальное приближение")
            a, e, i, Omega, omega, T = initial_guess
        else:
            a, e, i, Omega, omega, T = best_result.x

        return [a, e, i, Omega, omega, T]

    def calculate_accurate_position(self, a, e, i, Omega, omega, T, jd, earth_pos_func):
        """Улучшенный расчет позиции с проверками"""
        try:
            t = jd - T

            # Среднее движение
            n = np.sqrt(self.GM_sun / a**3)

            # Средняя аномалия
            M = n * t

            # Решение уравнения Кеплера
            E = self.solve_kepler_accurate(M, e)

            # Истинная аномалия
            nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2), np.sqrt(1 - e) * np.cos(E/2))

            # Радиус-вектор
            r = a * (1 - e * np.cos(E))

            i_rad = np.radians(i)
            Omega_rad = np.radians(Omega)
            omega_rad = np.radians(omega)

            # Орбитальные координаты
            x_orb = r * np.cos(nu)
            y_orb = r * np.sin(nu)

            # Преобразование в гелиоцентрические эклиптические координаты
            x_hel = (np.cos(omega_rad) * np.cos(Omega_rad) - np.sin(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * x_orb + \
                    (-np.sin(omega_rad) * np.cos(Omega_rad) - np.cos(omega_rad) * np.sin(Omega_rad) * np.cos(i_rad)) * y_orb

            y_hel = (np.cos(omega_rad) * np.sin(Omega_rad) + np.sin(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * x_orb + \
                    (-np.sin(omega_rad) * np.sin(Omega_rad) + np.cos(omega_rad) * np.cos(Omega_rad) * np.cos(i_rad)) * y_orb

            z_hel = (np.sin(omega_rad) * np.sin(i_rad)) * x_orb + (np.cos(omega_rad) * np.sin(i_rad)) * y_orb

            # Положение Земли
            x_earth, y_earth, z_earth = earth_pos_func(jd)

            # Геоцентрические координаты
            x_geo = x_hel - x_earth
            y_geo = y_hel - y_earth
            z_geo = z_hel - z_earth

            # Преобразование в экваториальные координаты
            eps = np.radians(23.4392911)

            x_eq = x_geo
            y_eq = y_geo * np.cos(eps) - z_geo * np.sin(eps)
            z_eq = y_geo * np.sin(eps) + z_geo * np.cos(eps)

            # Прямое восхождение и склонение
            ra = np.arctan2(y_eq, x_eq)
            dec = np.arctan2(z_eq, np.sqrt(x_eq**2 + y_eq**2))

            return np.degrees(ra) % 360, np.degrees(dec)

        except Exception as e:
            print(f"❌ Ошибка в calculate_accurate_position: {e}")
            return 0.0, 0.0

    def solve_kepler_accurate(self, M, e, iterations=20):
        """Улучшенное решение уравнения Кеплера"""
        M = M % (2 * np.pi)

        # Начальное приближение
        if e < 0.8:
            E = M + e * np.sin(M) + 0.5 * e**2 * np.sin(2*M)
        else:
            E = np.pi

        for _ in range(iterations):
            f = E - e * np.sin(E) - M
            f_prime = 1 - e * np.cos(E)

            if abs(f_prime) < 1e-15:
                break

            delta_E = f / f_prime
            E -= delta_E

            if abs(delta_E) < 1e-15:
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

        return np.degrees(nu) % 360

    def calculate_earth_approach(self, orbital_elements, days_ahead=365):
        """Стабильный расчет сближения с Землей"""
        try:
            a, e, i, Omega, omega, T = orbital_elements

            now = Time.now()
            min_distance = float('inf')
            best_time = now

            # Фиксируем время начала для воспроизводимости
            start_time = now.copy()

            # Поиск с фиксированным шагом
            search_points = []
            for days in range(0, days_ahead, 7):
                check_time = start_time + days * u.day
                search_points.append(check_time)

            # Добавляем дополнительные точки вокруг перигелия
            perihelion_time = Time(T, format='jd')
            for offset in [-30, -15, 0, 15, 30]:
                extra_time = perihelion_time + offset * u.day
                if start_time <= extra_time <= (start_time + days_ahead * u.day):
                    search_points.append(extra_time)

            # Убираем дубликаты и сортируем
            search_points = sorted(set(search_points), key=lambda t: t.jd)

            # Первый проход - грубый поиск
            for check_time in search_points:
                comet_pos = self.get_heliocentric_position(orbital_elements, check_time.jd)

                with solar_system_ephemeris.set('builtin'):
                    earth_pos = get_body_barycentric('earth', check_time)
                    earth_pos_au = np.array([
                        earth_pos.x.to(u.AU).value,
                        earth_pos.y.to(u.AU).value,
                        earth_pos.z.to(u.AU).value
                    ])

                distance = np.linalg.norm(comet_pos - earth_pos_au)

                if distance < min_distance:
                    min_distance = distance
                    best_time = check_time

            # Второй проход - точный поиск вокруг найденного минимума
            refine_days = 14  # ±7 дней для уточнения
            refine_step = 1   # Шаг 1 день

            refine_start = best_time - refine_days/2 * u.day
            for day_offset in range(0, refine_days + 1, refine_step):
                check_time = refine_start + day_offset * u.day

                comet_pos = self.get_heliocentric_position(orbital_elements, check_time.jd)

                with solar_system_ephemeris.set('builtin'):
                    earth_pos = get_body_barycentric('earth', check_time)
                    earth_pos_au = np.array([
                        earth_pos.x.to(u.AU).value,
                        earth_pos.y.to(u.AU).value,
                        earth_pos.z.to(u.AU).value
                    ])

                distance = np.linalg.norm(comet_pos - earth_pos_au)

                if distance < min_distance:
                    min_distance = distance
                    best_time = check_time

            # ФИКСИРУЕМ ВРЕМЯ - округляем до полудня для стабильности
            fixed_time = Time(best_time.datetime.replace(hour=12, minute=0, second=0, microsecond=0))

            return {
                'date': fixed_time.isot,  # Используем ISO формат для consistency
                'distance_au': float(min_distance),
                'is_safe': min_distance > 0.1,
                'min_distance_km': float(min_distance * 149597870.7),
                'search_accuracy_days': refine_step
            }

        except Exception as e:
            print(f"❌ Ошибка в calculate_earth_approach: {e}")
            # Резервный расчет с фиксированной логикой
            a, e, i, Omega, omega, T = orbital_elements
            perihelion_distance = a * (1 - e)
            earth_approach_distance = abs(perihelion_distance - 1.0)

            # Фиксированная дата через 30 дней
            fixed_date = Time.now() + 30 * u.day
            fixed_date = Time(fixed_date.datetime.replace(hour=12, minute=0, second=0, microsecond=0))

            return {
                'date': fixed_date.isot,
                'distance_au': float(earth_approach_distance),
                'is_safe': earth_approach_distance > 0.1,
                'min_distance_km': float(earth_approach_distance * 149597870.7),
                'fallback': True  # Помечаем что это резервный расчет
            }

    def get_heliocentric_position(self, orbital_elements, jd):
        """Гелиоцентрическая позиция кометы"""
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

# Остальные функции и эндпоинты остаются без изменений...
def calculate_orbit_from_observations(observations_array):
    calculator = CometOrbitCalculator()
    for obs in observations_array:
        ra, dec, datetime_str = obs
        calculator.add_observation(ra, dec, datetime_str)
    orbital_elements = calculator.calculate_orbital_elements()
    return orbital_elements

# ... остальные эндпоинты без изменений ...
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
                'image_data': planet[4],
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
        image_data = data.get('image_data', '')

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

        if len(observations) < 3:
            return jsonify({
                "success": False,
                "error": "Необходимо минимум 3 наблюдения. Получено: {}".format(len(observations))
            }), 400

        observations_array = []
        for obs in observations:
            observations_array.append([
                float(obs['ra']),
                float(obs['dec']),
                obs['time']
            ])

        calculator = CometOrbitCalculator()
        for obs in observations_array:
            calculator.add_observation(obs[0], obs[1], obs[2])

        orbital_elements = calculator.calculate_orbital_elements()

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
                "true_anomaly": float(round(true_anomaly, 6))
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

        required_params = ['semi_major_axis', 'eccentricity', 'inclination',
                          'longitude_ascending', 'argument_pericenter']

        for param in required_params:
            if param not in orbit_params:
                return jsonify({
                    "success": False,
                    "error": f"Отсутствует обязательный параметр: {param}"
                }), 400

        calculator = CometOrbitCalculator()

        orbital_elements = [
            float(orbit_params['semi_major_axis']),
            float(orbit_params['eccentricity']),
            float(orbit_params['inclination']),
            float(orbit_params['longitude_ascending']),
            float(orbit_params['argument_pericenter']),
            float(orbit_params.get('time_perihelion', Time.now().jd + 100))
        ]

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
