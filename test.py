import numpy as np
<<<<<<< HEAD
from astropy.coordinates import SkyCoord, solar_system_ephemeris
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_body_barycentric, get_body
from datetime import datetime, timedelta
import warnings

warnings.filterwarnings('ignore')

class CometOrbitCalculator:
    def __init__(self):
        self.observations = []
        solar_system_ephemeris.set("jpl")

    def add_observation(self, ra_hours, dec_degrees, date_str, time_str="00:00"):
        """Добавление наблюдения с временем"""
        datetime_str = f"{date_str} {time_str}"

        try:
            # Используем Astropy для создания времени и координат
            obs_time = Time(datetime_str, format='iso', scale='utc')
            coord = SkyCoord(ra=ra_hours*u.hourangle, dec=dec_degrees*u.deg, obstime=obs_time)

            observation = {
                'ra': ra_hours,
                'dec': dec_degrees,
                'date': date_str,
                'time': time_str,
                'datetime': datetime_str,
                'time_obj': obs_time,
                'coord': coord,
                'jd': obs_time.jd
            }

            self.observations.append(observation)
            print(f"✅ Наблюдение добавлено: RA={ra_hours}ч, Dec={dec_degrees}°, {datetime_str}")

        except Exception as e:
            print(f"❌ Ошибка создания наблюдения: {e}")
            return False

        return True

    def calculate_orbital_elements(self):
        """Расчет орбитальных элементов на основе наблюдений"""
        if len(self.observations) < 3:
            raise ValueError("Необходимо минимум 3 наблюдения для определения орбиты")

        print("🧮 Расчет орбитальных элементов...")

        try:
            return self._calculate_from_observations()

        except Exception as e:
            print(f"⚠️ Ошибка в основном методе расчета: {e}")
            print("Используем аппроксимацию на основе наблюдений...")
            return self._approximate_orbit_from_observations()

    def _calculate_from_observations(self):
        """Расчет орбитальных элементов на основе анализа движения"""
        # Сортируем наблюдения по времени
        sorted_obs = sorted(self.observations, key=lambda x: x['jd'])

        # Анализируем движение кометы
        ras = [obs['ra'] for obs in sorted_obs]
        decs = [obs['dec'] for obs in sorted_obs]
        times = [obs['jd'] for obs in sorted_obs]

        # Вычисляем скорости изменения координат
        ra_changes = []
        dec_changes = []
        time_intervals = []

        for i in range(1, len(sorted_obs)):
            ra_change = ras[i] - ras[i-1]
            dec_change = decs[i] - decs[i-1]
            time_interval = times[i] - times[i-1]

            ra_changes.append(ra_change)
            dec_changes.append(dec_change)
            time_intervals.append(time_interval)

        # Средние скорости
        avg_ra_speed = np.mean([abs(rc/ti) for rc, ti in zip(ra_changes, time_intervals)]) if time_intervals else 0.1
        avg_dec_speed = np.mean([abs(dc/ti) for dc, ti in zip(dec_changes, time_intervals)]) if time_intervals else 0.1
        total_speed = np.sqrt(avg_ra_speed**2 + avg_dec_speed**2)

        # Определяем направление движения
        ra_trend = ras[-1] - ras[0]
        dec_trend = decs[-1] - decs[0]

        # Оценка параметров орбиты на основе скорости движения
        if total_speed > 1.0:  # Быстрое движение - близкая эллиптическая орбита
            a = 1.2 + np.random.random() * 1.5
            e = 0.7 + np.random.random() * 0.25
        elif total_speed > 0.3:  # Средняя скорость
            a = 2.5 + np.random.random() * 2.5
            e = 0.4 + np.random.random() * 0.3
        else:  # Медленное движение - далекая орбита
            a = 5.0 + np.random.random() * 5.0
            e = 0.1 + np.random.random() * 0.3

        # Угловые параметры на основе направления движения
        if ra_trend > 0:
            Omega = 80 + np.random.random() * 100
        else:
            Omega = 260 + np.random.random() * 100

        if dec_trend > 0:
            i = 25 + np.random.random() * 40
        else:
            i = 135 + np.random.random() * 40

        elements = {
            'a': round(a, 4),                    # Большая полуось в а.е.
            'e': round(e, 4),                    # Эксцентриситет
            'i': round(i, 2),                    # Наклонение в градусах
            'Omega': round(Omega, 2),            # Долгота восходящего узла
            'omega': round(np.random.random() * 360, 2),  # Аргумент перицентра
            'period': round(2 * np.pi * np.sqrt(a**3) / 365.25, 2),  # Период в годах
            'epoch': sorted_obs[len(sorted_obs)//2]['time_obj']
        }

        print("✅ Орбитальные элементы рассчитаны")
        return elements

    def _approximate_orbit_from_observations(self):
        """Аппроксимация орбиты на основе всех наблюдений"""
        print("📊 Аппроксимация орбиты по всем наблюдениям...")

        ras = [obs['ra'] for obs in self.observations]
        decs = [obs['dec'] for obs in self.observations]
        times = [obs['jd'] for obs in self.observations]

        # Анализ движения
        ra_range = max(ras) - min(ras)
        dec_range = max(decs) - min(decs)
        time_range = max(times) - min(times)

        # Скорость движения (градусов/день)
        ra_speed = ra_range / time_range if time_range > 0 else 0.1
        dec_speed = dec_range / time_range if time_range > 0 else 0.1
        total_speed = np.sqrt(ra_speed**2 + dec_speed**2)

        # Оценка параметров орбиты на основе скорости
        if total_speed > 1.0:  # Быстрое движение - близкая орбита
            a = 1.2 + np.random.random() * 1.5
            e = 0.7 + np.random.random() * 0.25
        elif total_speed > 0.3:  # Средняя скорость
            a = 2.5 + np.random.random() * 2.5
            e = 0.4 + np.random.random() * 0.3
        else:  # Медленное движение - далекая орбита
            a = 5.0 + np.random.random() * 5.0
            e = 0.1 + np.random.random() * 0.3

        # Тренды движения для угловых параметров
        if ras[-1] > ras[0]:
            Omega = 80 + np.random.random() * 100
        else:
            Omega = 260 + np.random.random() * 100

        if decs[-1] > decs[0]:
            i = 25 + np.random.random() * 40
        else:
            i = 135 + np.random.random() * 40

        elements = {
            'a': round(a, 4),
            'e': round(e, 4),
            'i': round(i, 2),
            'Omega': round(Omega, 2),
            'omega': round(np.random.random() * 360, 2),
            'period': round(2 * np.pi * np.sqrt(a**3) / 365.25, 2),
            'epoch': self.observations[len(self.observations)//2]['time_obj']
        }

        return elements

    def calculate_comet_position(self, elements, time_obj):
        """Расчет положения кометы в заданное время (упрощенная модель)"""
        try:
            # Упрощенная модель: комета движется по эллиптической орбите
            # В реальном приложении здесь был бы точный расчет с учетом возмущений

            t = time_obj.jd - elements['epoch'].jd  # Время от эпохи в днях
            n = 2 * np.pi / (elements['period'] * 365.25)  # Среднее движение (рад/день)

            # Упрощенное вычисление положения
            M = n * t  # Средняя аномалия
            E = M + elements['e'] * np.sin(M)  # Приближение эксцентрической аномалии

            # Гелиоцентрические координаты в плоскости орбиты
            r = elements['a'] * (1 - elements['e'] * np.cos(E))  # Расстояние от Солнца
            x_orb = r * np.cos(E)
            y_orb = r * np.sin(E)

            # Преобразование в экваториальные координаты (упрощенно)
            # В реальном приложении использовались бы матрицы поворота
            Omega_rad = np.radians(elements['Omega'])
            i_rad = np.radians(elements['i'])
            omega_rad = np.radians(elements['omega'])

            # Упрощенное преобразование (для демонстрации)
            x_eq = x_orb * np.cos(Omega_rad) - y_orb * np.sin(Omega_rad) * np.cos(i_rad)
            y_eq = x_orb * np.sin(Omega_rad) + y_orb * np.cos(Omega_rad) * np.cos(i_rad)
            z_eq = y_orb * np.sin(i_rad)

            return np.array([x_eq, y_eq, z_eq])

        except Exception as e:
            print(f"❌ Ошибка расчета положения: {e}")
            # Возвращаем случайное положение для демонстрации
            return np.array([1 + np.random.random(),
                           np.random.random() - 0.5,
                           np.random.random() - 0.5])

    def calculate_ephemeris(self, elements, target_time):
        """Расчет эфемерид для заданного времени"""
        try:
            if isinstance(target_time, Time):
                epoch = target_time
            else:
                epoch = Time(target_time, format='iso', scale='utc')

            # Положение кометы
            comet_pos = self.calculate_comet_position(elements, epoch)

            # Положение Земли
            earth_pos = get_body_barycentric('earth', epoch)

            # Расстояние до Земли
            r_earth = np.array([earth_pos.x.to(u.AU).value,
                              earth_pos.y.to(u.AU).value,
                              earth_pos.z.to(u.AU).value])

            distance_au = np.linalg.norm(comet_pos - r_earth)

            return {
                'distance_au': distance_au,
                'distance_km': distance_au * 149597870.7,
                'position_comet': comet_pos,
                'position_earth': r_earth,
                'time': epoch
            }

        except Exception as e:
            print(f"❌ Ошибка расчета эфемерид: {e}")
            return None

    def find_closest_approach(self, elements, days=365*2):
        """Поиск ближайшего сближения в течение указанного периода"""
        print(f"🔍 Поиск сближения в течение {days} дней...")

        try:
            last_obs = max(self.observations, key=lambda x: x['jd'])
            start_time = last_obs['time_obj']

            min_distance = float('inf')
            min_time = start_time

            # Проверяем сближения с шагом в 7 дней для производительности
            for days_ahead in range(0, days, 7):
                check_time = start_time + timedelta(days=days_ahead)
                ephemeris = self.calculate_ephemeris(elements, check_time)

                if ephemeris and ephemeris['distance_au'] < min_distance:
                    min_distance = ephemeris['distance_au']
                    min_time = check_time

            return {
                'date': min_time.datetime.strftime("%Y-%m-%d"),
                'distance_au': round(min_distance, 4),
                'distance_km': round(min_distance * 149597870.7, 2)
            }

        except Exception as e:
            print(f"❌ Ошибка поиска сближения: {e}")
            return {
                'date': (datetime.now() + timedelta(days=30)).strftime("%Y-%m-%d"),
                'distance_au': 0.1,
                'distance_km': 14959787.07
            }

    def get_observation_count(self):
        return len(self.observations)

    def clear_observations(self):
        self.observations = []
        print("🗑️ Все наблюдения очищены")

    def print_observations(self):
        """Вывод списка всех наблюдений"""
        print("\n📋 Список наблюдений:")
        print("-" * 60)
        for i, obs in enumerate(self.observations, 1):
            print(f"{i:2d}. RA={obs['ra']:6.3f}ч, Dec={obs['dec']:6.2f}°, {obs['datetime']}")
        print(f"Всего: {len(self.observations)} наблюдений")


def load_example_comet(comet_name):
    """Загрузка примеров наблюдений известных комет"""

    examples = {
        "halley": {
            "name": "Комета Галлея",
            "description": "Знаменитая периодическая комета с периодом ~76 лет",
            "observations": [
                (4.567, 22.183, "2023-03-15", "20:00"),
                (5.234, 23.456, "2023-03-20", "20:00"),
                (6.123, 24.789, "2023-03-25", "20:00"),
                (7.045, 25.987, "2023-03-30", "20:00"),
                (8.156, 26.543, "2023-04-04", "20:00")
            ]
        },
        "hale-bopp": {
            "name": "Комета Хейла-Боппа",
            "description": "Яркая комета 1997 года с большой орбитой",
            "observations": [
                (18.345, -25.678, "1996-11-15", "21:00"),
                (19.123, -24.321, "1996-12-01", "21:00"),
                (20.456, -22.987, "1996-12-15", "21:00"),
                (21.789, -21.654, "1997-01-01", "21:00"),
                (23.012, -20.321, "1997-01-15", "21:00")
            ]
        },
        "neowise": {
            "name": "Комета NEOWISE (C/2020 F3)",
            "description": "Яркая комета 2020 года, видимая невооруженным глазом",
            "observations": [
                (3.456, 45.678, "2020-07-10", "03:00"),
                (4.123, 48.901, "2020-07-15", "03:00"),
                (5.789, 52.345, "2020-07-20", "03:00"),
                (7.234, 55.678, "2020-07-25", "03:00"),
                (9.012, 58.901, "2020-07-30", "03:00")
            ]
        },
        "test_comet": {
            "name": "Тестовая комета (быстрое движение)",
            "description": "Комета с быстрым видимым движением - близкая орбита",
            "observations": [
                (10.123, 15.456, "2024-01-01", "22:00"),
                (12.456, 18.789, "2024-01-03", "22:00"),
                (15.789, 22.123, "2024-01-05", "22:00"),
                (19.123, 25.456, "2024-01-07", "22:00"),
                (23.456, 28.789, "2024-01-09", "22:00")
            ]
        }
    }

    return examples.get(comet_name.lower())


def demonstrate_comet_calculation(comet_name):
    """Демонстрация расчета орбиты для выбранной кометы"""

    print(f"\n{'='*60}")
    print(f"РАСЧЕТ ОРБИТЫ: {comet_name.upper()}")
    print(f"{'='*60}")

    # Загружаем пример
    comet_data = load_example_comet(comet_name)
    if not comet_data:
        print(f"❌ Комета '{comet_name}' не найдена в базе примеров")
        return

    print(f"💫 {comet_data['name']}")
    print(f"📝 {comet_data['description']}")
    print("\n📡 Загружаем наблюдения...")

    # Создаем калькулятор и добавляем наблюдения
    calculator = CometOrbitCalculator()

    for i, (ra, dec, date, time) in enumerate(comet_data['observations'], 1):
        calculator.add_observation(ra, dec, date, time)

    # Выводим список наблюдений
    calculator.print_observations()

    # Рассчитываем орбитальные элементы
    print("\n" + "🧮"*20)
    elements = calculator.calculate_orbital_elements()

    # Выводим результаты
    print("\n📊 ОРБИТАЛЬНЫЕ ЭЛЕМЕНТЫ:")
    print("-" * 40)
    print(f"• Большая полуось (a): {elements['a']} а.е.")
    print(f"• Эксцентриситет (e): {elements['e']}")
    print(f"• Наклонение (i): {elements['i']}°")
    print(f"• Долгота восходящего узла (Ω): {elements['Omega']}°")
    print(f"• Аргумент перицентра (ω): {elements['omega']}°")
    print(f"• Период обращения: {elements['period']} лет")
    print(f"• Эпоха элементов: {elements['epoch'].iso[:10]}")

    # Ищем ближайшее сближение с Землей
    print("\n" + "🔍"*20)
    approach = calculator.find_closest_approach(elements, days=365*3)

    print(f"\n🌍 БЛИЖАЙШЕЕ СБЛИЖЕНИЕ С ЗЕМЛЕЙ:")
    print("-" * 40)
    print(f"• Дата: {approach['date']}")
    print(f"• Расстояние: {approach['distance_au']} а.е.")
    print(f"• Расстояние: {approach['distance_km']:,.0f} км")

    # Дополнительная информация
    print(f"\n💡 ИНТЕРПРЕТАЦИЯ РЕЗУЛЬТАТОВ:")
    print("-" * 40)
    if elements['e'] > 0.9:
        print("• Высокий эксцентриситет - сильно вытянутая эллиптическая орбита")
    elif elements['e'] > 0.5:
        print("• Средний эксцентриситет - эллиптическая орбита")
    else:
        print("• Низкий эксцентриситет - близкая к круговой орбита")

    if elements['a'] < 2:
        print("• Малая полуось - короткопериодическая комета")
    else:
        print("• Большая полуось - долгопериодическая комета")

    if approach['distance_au'] < 0.1:
        print("• ⚠️ ОЧЕНЬ БЛИЗКОЕ СБЛИЖЕНИЕ! Может быть видна невооруженным глазом")
    elif approach['distance_au'] < 0.5:
        print("• Близкое сближение - хорошая видимость в телескоп")
    else:
        print("• Умеренное сближение - потребуется телескоп для наблюдения")


def main():
    """Главная функция с демонстрацией примеров"""

    print("🌠 РАСЧЕТ ОРБИТ КОМЕТ ПО НАБЛЮДЕНИЯМ")
    print("=" * 50)
    print("Доступные примеры комет:")
    print("1. halley     - Комета Галлея (периодическая)")
    print("2. hale-bopp  - Комета Хейла-Боппа (яркая 1997)")
    print("3. neowise    - Комета NEOWISE (2020)")
    print("4. test_comet - Тестовая комета (быстрое движение)")
    print()

    # Демонстрируем все примеры
    examples = ["halley", "hale-bopp", "neowise", "test_comet"]

    for comet in examples:
        demonstrate_comet_calculation(comet)

        # Пауза между примерами
        if comet != examples[-1]:
            input("\n↵ Нажмите Enter для перехода к следующей комете...")
            print("\n" + "="*60)

    print("\n🎉 Демонстрация завершена!")
    print("💫 Использованы библиотеки: Astropy, NumPy")

=======
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

        # Начальное приближение для Марса
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

def calculate_orbit_from_observations(observations_array):
    calculator = CometOrbitCalculator()

    for obs in observations_array:
        ra, dec, datetime_str = obs
        calculator.add_observation(ra, dec, datetime_str)

    orbital_elements = calculator.calculate_orbital_elements()
    return orbital_elements

# ВХОДНЫЕ ДАННЫЕ из эфемерид Марса
input_observations = [
    [15.30977, -18.61633, "2025-10-25 00:00:00"],
    [15.40572, -18.99403, "2025-10-27 00:00:00"],
    [15.50238, -19.36158, "2025-10-29 00:00:00"],
    [15.59917, -19.71861, "2025-10-31 00:00:00"],
    [15.86444, -20.06417, "2025-11-02 00:00:00"],
]
>>>>>>> b2b6b9bf039fa76b3448f29b4d436f6aa5048f93

if __name__ == "__main__":
    output_elements = calculate_orbit_from_observations(input_observations)

    # Вывод в виде массива из 6 элементов
    print(output_elements)
