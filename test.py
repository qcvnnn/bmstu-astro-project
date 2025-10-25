import numpy as np
from scipy.optimize import fsolve
import math
from datetime import datetime, timedelta

class CometOrbitCalculator:
    def __init__(self):
        # Астрономические константы
        self.k = 0.01720209895  # Gaussian gravitational constant
        self.au_km = 149597870.7  # 1 а.е. в км

        # Хранилище наблюдений
        self.observations = []

    def add_observation(self, ra_hours, dec_degrees, date_str, time_str="00:00"):
        """Добавление наблюдения с временем"""
        # Объединяем дату и время
        datetime_str = f"{date_str} {time_str}"

        observation = {
            'ra': ra_hours,
            'dec': dec_degrees,
            'date': date_str,
            'time': time_str,
            'datetime': datetime_str,
            'jd': self.datetime_to_jd(datetime_str)
        }

        self.observations.append(observation)
        return len(self.observations)

    def datetime_to_jd(self, datetime_str):
        """Преобразование даты и времени в Юлианскую дату"""
        try:
            # Парсим строку даты-времени
            dt = datetime.strptime(datetime_str, "%Y-%m-%d %H:%M")
        except ValueError:
            try:
                dt = datetime.strptime(datetime_str, "%Y-%m-%d")
            except ValueError:
                # Если не удается распарсить, используем только дату
                dt = datetime.strptime(datetime_str.split()[0], "%Y-%m-%d")

        # Формула для преобразования в Юлианскую дату
        a = (14 - dt.month) // 12
        y = dt.year + 4800 - a
        m = dt.month + 12 * a - 3

        # Целая часть Юлианской даты
        jd = dt.day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045

        # Добавляем дробную часть (время суток)
        time_fraction = (dt.hour + dt.minute / 60.0) / 24.0
        jd += time_fraction

        return jd

    def jd_to_datetime(self, jd):
        """Преобразование Юлианской даты в datetime"""
        jd = jd + 0.5
        Z = int(jd)
        F = jd - Z

        if Z < 2299161:
            A = Z
        else:
            alpha = int((Z - 1867216.25) / 36524.25)
            A = Z + 1 + alpha - int(alpha / 4)

        B = A + 1524
        C = int((B - 122.1) / 365.25)
        D = int(365.25 * C)
        E = int((B - D) / 30.6001)

        day = int(B - D - int(30.6001 * E) + F)
        month = E - 1 if E < 14 else E - 13
        year = C - 4716 if month > 2 else C - 4715

        # Вычисляем время
        time_fraction = F - int(F)
        hours = int(time_fraction * 24)
        minutes = int((time_fraction * 24 - hours) * 60)

        return datetime(year, month, day, hours, minutes)

    def clear_observations(self):
        """Очистка всех наблюдений"""
        self.observations = []

    def get_observation_count(self):
        """Получить количество наблюдений"""
        return len(self.observations)

    def ra_dec_to_equatorial(self, ra_hours, dec_degrees):
        """Преобразование прямого восхождения и склонения в экваториальные координаты"""
        # Преобразование часов в градусы
        ra_degrees = ra_hours * 15.0  # 1 час = 15 градусов

        ra_rad = np.radians(ra_degrees)
        dec_rad = np.radians(dec_degrees)

        # Единичный вектор направления
        x = np.cos(dec_rad) * np.cos(ra_rad)
        y = np.cos(dec_rad) * np.sin(ra_rad)
        z = np.sin(dec_rad)

        return np.array([x, y, z])

    def calculate_orbital_elements(self):
        """Расчет орбитальных элементов на основе всех наблюдений"""
        if len(self.observations) < 3:
            raise ValueError("Необходимо как минимум 3 наблюдения для расчета орбиты")

        # Сортируем наблюдения по времени
        sorted_obs = sorted(self.observations, key=lambda x: x['jd'])

        # Используем метод Гаусса с тремя наблюдениями
        return self.gauss_method(sorted_obs)

    def gauss_method(self, observations):
        """Метод Гаусса для определения орбитальных элементов"""
        # Выбираем три наблюдения (первое, среднее и последнее)
        n = len(observations)
        idx1, idx2, idx3 = 0, n//2, n-1

        obs1, obs2, obs3 = observations[idx1], observations[idx2], observations[idx3]

        ra1, ra2, ra3 = obs1['ra'], obs2['ra'], obs3['ra']
        dec1, dec2, dec3 = obs1['dec'], obs2['dec'], obs3['dec']
        jd1, jd2, jd3 = obs1['jd'], obs2['jd'], obs3['jd']

        # Вычисляем тау-интервалы
        tau1 = self.k * (jd3 - jd2)
        tau3 = self.k * (jd2 - jd1)
        tau = self.k * (jd3 - jd1)

        # Единичные векторы направлений
        rho1 = self.ra_dec_to_equatorial(ra1, dec1)
        rho2 = self.ra_dec_to_equatorial(ra2, dec2)
        rho3 = self.ra_dec_to_equatorial(ra3, dec3)

        # Вычисляем вспомогательные величины
        D0 = np.dot(rho1, np.cross(rho2, rho3))
        D11 = np.dot(np.cross(rho1, rho2), rho3)
        D21 = np.dot(np.cross(rho1, rho3), rho2)
        D31 = np.dot(np.cross(rho2, rho3), rho1)

        # Вычисляем геоцентрические расстояния
        A1 = tau3 / tau
        B1 = A1 * (tau**2 - tau3**2) / 6.0
        A3 = -tau1 / tau
        B3 = A3 * (tau**2 - tau1**2) / 6.0

        # Решаем систему уравнений для расстояний
        def equations(vars):
            r2, rho1_mag, rho3_mag = vars

            eq1 = A1 * (rho1_mag * D21 / D0 + 1/r2**3 * B1 * D21 / D0) + \
                  rho3_mag * D31 / D0 + 1/r2**3 * B3 * D31 / D0 - r2

            eq2 = rho1_mag - (D11 / D0) * (1 + B1 / r2**3)
            eq3 = rho3_mag - (D31 / D0) * (1 + B3 / r2**3)

            return [eq1, eq2, eq3]

        # Начальное приближение
        r2_guess = 2.0  # а.е.
        rho1_guess = 0.1
        rho3_guess = 0.1

        solution = fsolve(equations, [r2_guess, rho1_guess, rho3_guess])
        r2, rho1_mag, rho3_mag = solution

        # Вычисляем гелиоцентрические положения (упрощенно)
        R1 = np.array([0, 0, 0])  # Положение Земли (упрощенно)
        R2 = np.array([0, 0, 0])
        R3 = np.array([0, 0, 0])

        r1_vec = rho1_mag * rho1 - R1
        r2_vec = r2 * rho2 - R2
        r3_vec = rho3_mag * rho3 - R3

        # Вычисляем скорость (упрощенно)
        v2_vec = (r3_vec - r1_vec) / (tau1 + tau3)

        # Вычисляем орбитальные элементы
        return self.vectors_to_orbital_elements(r2_vec, v2_vec)

    def vectors_to_orbital_elements(self, r, v):
        """Преобразование векторов положения и скорости в орбитальные элементы"""
        mu = self.k**2  # Гравитационный параметр

        # Удельный момент импульса
        h = np.cross(r, v)
        h_mag = np.linalg.norm(h)

        # Вектор эксцентриситета
        r_mag = np.linalg.norm(r)
        v_mag = np.linalg.norm(v)

        e_vec = ((v_mag**2 - mu/r_mag) * r - np.dot(r, v) * v) / mu
        e = np.linalg.norm(e_vec)

        # Большая полуось
        energy = v_mag**2 / 2 - mu / r_mag
        if abs(energy) < 1e-10:
            a = float('inf')  # Параболическая орбита
        else:
            a = -mu / (2 * energy) if energy < 0 else mu / (2 * energy)

        # Наклонение
        i = np.degrees(np.arccos(h[2] / h_mag))

        # Долгота восходящего узла
        node_vec = np.cross([0, 0, 1], h)
        node_mag = np.linalg.norm(node_vec)
        if node_mag > 1e-10:
            Omega = np.degrees(np.arctan2(node_vec[0], node_vec[1]))
            if Omega < 0:
                Omega += 360
        else:
            Omega = 0

        # Аргумент перицентра
        if e > 1e-10 and node_mag > 1e-10:
            n = node_vec / node_mag
            cos_omega = np.dot(n, e_vec) / (node_mag * e)
            cos_omega = np.clip(cos_omega, -1, 1)
            omega = np.degrees(np.arccos(cos_omega))
            if e_vec[2] < 0:
                omega = 360 - omega
        else:
            omega = 0

        # Истинная аномалия
        if e > 1e-10:
            cos_nu = np.dot(e_vec, r) / (e * r_mag)
            cos_nu = np.clip(cos_nu, -1, 1)
            nu = np.degrees(np.arccos(cos_nu))
            if np.dot(r, v) < 0:
                nu = 360 - nu
        else:
            nu = np.degrees(np.arctan2(r[1], r[0]))

        # Период (для эллиптических орбит)
        if a > 0 and not math.isinf(a):
            period = 2 * np.pi * np.sqrt(a**3 / mu) / 365.25  # в годах
        else:
            period = float('inf')

        return {
            'a': a,
            'e': e,
            'i': i,
            'Omega': Omega,
            'omega': omega,
            'nu': nu,
            'period': period
        }

    def calculate_ephemeris(self, elements, jd):
        """Расчет эфемерид для заданной даты"""
        try:
            # Упрощенный расчет положения
            a, e, nu = elements['a'], elements['e'], elements['nu']

            # Расчет расстояния
            r = a * (1 - e**2) / (1 + e * np.cos(np.radians(nu)))

            return {
                'distance_au': r,
                'distance_km': r * self.au_km
            }
        except:
            return {'distance_au': 0.1, 'distance_km': 14959787.0}

    def find_closest_approach(self, elements, start_date=None, days=365):
        """Поиск ближайшего сближения с Землей"""
        if start_date is None:
            # Используем дату последнего наблюдения
            last_obs = max(self.observations, key=lambda x: x['jd'])
            start_jd = last_obs['jd']
        else:
            start_jd = self.datetime_to_jd(start_date + " 00:00")

        min_distance = float('inf')
        min_date = start_date
        min_jd = start_jd

        for day in range(days):
            current_jd = start_jd + day
            ephemeris = self.calculate_ephemeris(elements, current_jd)

            if ephemeris['distance_au'] < min_distance:
                min_distance = ephemeris['distance_au']
                min_jd = current_jd

        min_datetime = self.jd_to_datetime(min_jd)
        return min_datetime.strftime("%Y-%m-%d %H:%M"), min_distance


def main():
    """Пример использования с временем"""
    print("Введите 5 наблюдений (формат: ЧЧ.Ч ГГ.Г ГГГГ-ММ-ДД ЧЧ:ММ)")
    print("Если время не указано, используется 00:00")
    print()

    calculator = CometOrbitCalculator()

    for i in range(5):
        print(f"--- Наблюдение {i+1} ---")
        try:
            ra = float(input("Прямое восхождение (часы): "))
            dec = float(input("Склонение (градусы): "))
            date = input("Дата (ГГГГ-ММ-ДД): ")
            time = input("Время (ЧЧ:ММ) [опционально]: ").strip()

            if not time:
                time = "00:00"

            calculator.add_observation(ra, dec, date, time)
            print("📌 Добавлено\n")

        except ValueError as e:
            print(f"❌ Ошибка ввода: {e}")
            return
        except Exception as e:
            print(f"❌ Ошибка: {e}")
            return

    print("📌 Расчет...")

    try:
        # Вычисляем орбитальные элементы
        elements = calculator.calculate_orbital_elements()

        # Находим ближайшее сближение
        approach_datetime, approach_dist = calculator.find_closest_approach(elements)

        print("\n📌 Результаты:")
        print(f"a: {elements['a']:.3f} a.e.")
        print(f"e: {elements['e']:.4f}")
        print(f"i: {elements['i']:.2f}°")
        print(f"Ω: {elements['Omega']:.1f}°")
        print(f"ω: {elements['omega']:.1f}°")
        print(f"Период: {elements['period']:.2f} лет")
        print(f"Сближение: {approach_datetime}")
        print(f"Расстояние: {approach_dist:.3f} a.e.")

    except Exception as e:
        print(f"❌ Ошибка расчета: {e}")
        import traceback
        traceback.print_exc()




def validate_datetime(date_str, time_str):
    """Проверка корректности даты и времени"""
    try:
        datetime_str = f"{date_str} {time_str}"
        datetime.strptime(datetime_str, "%Y-%m-%d %H:%M")
        return True, ""
    except ValueError as e:
        return False, f"Неверный формат: {e}"

def manual_input_5chars():
    """Ручной ввод с массивами по 5 символов с поддержкой времени"""
    print("🎯 РУЧНОЙ ВВОД НАБЛЮДЕНИЙ")
    print("=" * 40)
    print("Формат: ЧЧ.Ч ГГ.Г ГГГГ-ММ-ДД ЧЧ:ММ")
    print("Если время не указано, используется 00:00")

    calculator = CometOrbitCalculator()

    for i in range(5):
        print(f"\n--- Наблюдение {i+1} ---")

        while True:
            ra = input("Прямое восхождение (часы): ").strip()[:5]
            dec = input("Склонение (градусы): ").strip()[:5]
            date_str = input("Дата (ГГГГ-ММ-ДД): ").strip()[:10]
            time_str = input("Время (ЧЧ:ММ): ").strip()[:5]

            # Если время не указано, используем 00:00
            if not time_str:
                time_str = "00:00"
            # Если указано без минут, добавляем :00
            elif ":" not in time_str and len(time_str) <= 2:
                time_str = f"{time_str}:00"

            try:
                # Проверяем числа
                ra_val = float(ra)
                dec_val = float(dec)

                # Проверяем дату и время
                valid, error_msg = validate_datetime(date_str, time_str)
                if not valid:
                    print(f"❌ Ошибка в дате/времени: {error_msg}")
                    print("Попробуйте снова\n")
                    continue

                # Если все проверки пройдены
                calculator.add_observation(ra_val, dec_val, date_str, time_str)
                print("✅ Добавлено")
                break

            except ValueError as e:
                print(f"❌ Ошибка в числах: {e}")
                print("Попробуйте снова\n")

    return calculator, 5

def main():
    calculator, count = manual_input_5chars()

    if count < 3:
        print(f"\nНедостаточно данных: {count}")
        return

    print("\n🧮 Расчет...")

    try:
        elements = calculator.calculate_orbital_elements()

        # Находим последнее наблюдение для старта поиска сближения
        last_observation = max(calculator.observations, key=lambda x: x['jd'])
        approach_date, approach_dist = calculator.find_closest_approach(elements)

        print("\n📊 Результаты:")
        print(f"Большая полуось: {elements['a']:.3f} а.е.")
        print(f"Эксцентриситет: {elements['e']:.4f}")
        print(f"Наклонение: {elements['i']:.2f}°")
        print(f"Долгота узла: {elements['Omega']:.1f}°")
        print(f"Аргумент перицентра: {elements['omega']:.1f}°")
        print(f"Период: {elements['period']:.2f} лет")
        print(f"Сближение: {approach_date}")
        print(f"Расстояние: {approach_dist:.3f} а.е.")

        # Дополнительная информация о наблюдениях
        print(f"\n📅 Использовано наблюдений: {calculator.get_observation_count()}")
        print("Временные метки наблюдений:")
        for i, obs in enumerate(calculator.observations, 1):
            print(f"  {i}. {obs['datetime']} (JD: {obs['jd']:.4f})")

    except Exception as e:
        print(f"Ошибка расчета: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
