"""
Ручной ввод наблюдений - массивы по 5 символов с поддержкой времени
"""

from orbit_calculator import CometOrbitCalculator
from datetime import datetime

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
