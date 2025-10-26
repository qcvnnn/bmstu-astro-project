let pointCounter = 5;
let currentOrbitData = null;

// Функции для работы с изображениями
function initImageUpload() {
    const imageUpload = document.getElementById('imageUpload');
    const imagePreview = document.getElementById('imagePreview');

    imageUpload.addEventListener('change', handleImageUpload);

    // Добавляем обработчик перетаскивания
    imagePreview.addEventListener('dragover', function(e) {
        e.preventDefault();
        this.style.borderColor = '#3498db';
        this.style.background = '#f8f9fa';
    });

    imagePreview.addEventListener('dragleave', function(e) {
        e.preventDefault();
        this.style.borderColor = '#bdc3c7';
        this.style.background = '#ffffff';
    });

    imagePreview.addEventListener('drop', function(e) {
        e.preventDefault();
        this.style.borderColor = '#bdc3c7';
        this.style.background = '#ffffff';

        const files = e.dataTransfer.files;
        if (files.length > 0) {
            imageUpload.files = files;
            handleImageUpload({ target: imageUpload });
        }
    });
}

function handleImageUpload(event) {
    const file = event.target.files[0];
    if (!file) return;

    // Проверяем тип файла
    if (!file.type.match('image.*')) {
        alert('Пожалуйста, выберите файл изображения (JPEG, PNG, GIF и т.д.)');
        return;
    }

    // Проверяем размер файла (максимум 5MB)
    if (file.size > 5 * 1024 * 1024) {
        alert('Размер файла не должен превышать 5MB');
        return;
    }

    const reader = new FileReader();

    reader.onload = function(e) {
        const imagePreview = document.getElementById('imagePreview');
        imagePreview.innerHTML = '';

        const img = document.createElement('img');
        img.src = e.target.result;
        img.alt = 'Изображение кометы';

        imagePreview.appendChild(img);

        // Сохраняем изображение в localStorage
        localStorage.setItem('cometImage', e.target.result);

        showNotification('✅ Изображение успешно загружено!', 'success');
    };

    reader.onerror = function() {
        showNotification('❌ Ошибка при чтении файла', 'error');
    };

    reader.readAsDataURL(file);
}

function removeImage() {
    const imagePreview = document.getElementById('imagePreview');
    const imageUpload = document.getElementById('imageUpload');

    imagePreview.innerHTML = `
        <div class="placeholder-content">
            <div class="placeholder-icon"></div>
            <p>Перетащите сюда изображение<br>или нажмите кнопку ниже</p>
        </div>
    `;
    imageUpload.value = '';

    // Удаляем из localStorage
    localStorage.removeItem('cometImage');

    showNotification('🗑 Изображение удалено', 'info');
}

function restoreImage() {
    const savedImage = localStorage.getItem('cometImage');
    if (savedImage) {
        const imagePreview = document.getElementById('imagePreview');
        imagePreview.innerHTML = '';

        const img = document.createElement('img');
        img.src = savedImage;
        img.alt = 'Изображение кометы';

        imagePreview.appendChild(img);
    }
}

function showNotification(message, type) {
    // Создаем уведомление
    const notification = document.createElement('div');
    notification.style.cssText = `
        position: fixed;
        top: 20px;
        right: 20px;
        padding: 15px 20px;
        background: ${type === 'success' ? '#27ae60' : type === 'error' ? '#e74c3c' : '#3498db'};
        color: white;
        border-radius: 10px;
        box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        z-index: 1000;
        font-weight: 600;
        transform: translateX(100%);
        transition: transform 0.3s ease;
    `;
    notification.textContent = message;

    document.body.appendChild(notification);

    // Анимация появления
    setTimeout(() => {
        notification.style.transform = 'translateX(0)';
    }, 100);
        // Автоматическое скрытие
        setTimeout(() => {
          notification.style.transform = 'translateX(100%)';
          setTimeout(() => {
              document.body.removeChild(notification);
          }, 300);
      }, 3000);
  }

// ФУНКЦИИ ДЛЯ РАБОТЫ С БАЗОЙ ДАННЫХ ПЛАНЕТ

async function loadPlanets() {
    try {
        const response = await fetch('http://127.0.0.1:5001/api/planets');
        const result = await response.json();

        if (result.success) {
            displayPlanets(result.planets);
        } else {
            alert('Ошибка загрузки планет: ' + result.error);
        }
    } catch (error) {
        alert('Ошибка соединения: ' + error.message);
    }
}

function displayPlanets(planets) {
  const planetsList = document.getElementById('planets-list');
  planetsList.innerHTML = '';

  if (planets.length === 0) {
      planetsList.innerHTML = '<p style="text-align: center; color: #666; padding: 20px;">Нет сохраненных планет</p>';
      return;
  }

  planets.forEach(planet => {
      const planetElement = document.createElement('div');
      planetElement.className = 'planet-card';

      // ДОБАВЛЯЕМ ПРЕВЬЮ ИЗОБРАЖЕНИЯ
      const imagePreview = planet.image_data ?
          `<div class="planet-image-preview">
              <img src="${planet.image_data}" alt="${planet.name}" onclick="showFullImage('${planet.image_data}')">
          </div>` :
          '<div class="planet-no-image">📷 Нет изображения</div>';

      planetElement.innerHTML = `
          <div class="planet-header">
              <h3>${planet.name}</h3>
              <button class="delete-btn" onclick="deletePlanet(${planet.id})">🗑️ Удалить</button>
          </div>
          ${imagePreview}
          <div class="planet-info">
              <p><strong>Наблюдения:</strong> ${planet.observations.length} точек</p>
              <p><strong>Большая полуось:</strong> ${planet.orbital_elements.semi_major_axis} а.е.</p>
              <p><strong>Эксцентриситет:</strong> ${planet.orbital_elements.eccentricity}</p>
              <p><strong>Создана:</strong> ${new Date(planet.created_at).toLocaleString()}</p>
          </div>
          <button class="load-btn" onclick="loadPlanetData(${planet.id})">📊 Загрузить данные</button>
      `;
      planetsList.appendChild(planetElement);
  });
}

// ФУНКЦИЯ ДЛЯ ПОКАЗА ИЗОБРАЖЕНИЯ В ПОЛНОМ РАЗМЕРЕ
function showFullImage(imageData) {
  const modal = document.createElement('div');
  modal.style.cssText = `
      position: fixed;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background: rgba(0,0,0,0.8);
      display: flex;
      justify-content: center;
      align-items: center;
      z-index: 1000;
      cursor: pointer;
  `;

  const img = document.createElement('img');
  img.src = imageData;
  img.style.cssText = `
      max-width: 90%;
      max-height: 90%;
      object-fit: contain;
      border-radius: 10px;
  `;

  modal.appendChild(img);
  modal.onclick = () => document.body.removeChild(modal);
  document.body.appendChild(modal);
}

async function deletePlanet(planetId) {
    if (!confirm('Удалить эту планету?')) return;

    try {
        const response = await fetch(`http://127.0.0.1:5001/api/planets/${planetId}`, {
            method: 'DELETE'
        });
        const result = await response.json();

        if (result.success) {
            alert('Планета удалена');
            loadPlanets();
        } else {
            alert('Ошибка удаления: ' + result.error);
        }
    } catch (error) {
        alert('Ошибка соединения: ' + error.message);
    }
}

async function loadPlanetData(planetId) {
  try {
      const response = await fetch('http://127.0.0.1:5001/api/planets');
      const result = await response.json();

      if (result.success) {
          const planet = result.planets.find(p => p.id === planetId);
          if (planet) {
              // Заполняем поля наблюдениями
              fillObservations(planet.observations);
              // Заполняем результаты орбиты
              fillOrbitResults(planet.orbital_elements);
              // ЗАГРУЖАЕМ ИЗОБРАЖЕНИЕ
              if (planet.image_data) {
                  loadPlanetImage(planet.image_data);
              } else {
                  removeImage(); // Очищаем если нет изображения
              }
              alert(`✅ Данные планеты "${planet.name}" загружены`);
          }
      }
  } catch (error) {
      alert('Ошибка загрузки данных: ' + error.message);
  }
}

function loadPlanetImage(imageData) {
  const imagePreview = document.getElementById('imagePreview');
  if (imageData && imagePreview) {
      imagePreview.innerHTML = '';
      const img = document.createElement('img');
      img.src = imageData;
      img.alt = 'Изображение кометы';
      imagePreview.appendChild(img);

      // Сохраняем в localStorage для текущей сессии
      localStorage.setItem('cometImage', imageData);
  }
}

function fillObservations(observations) {
    // Очищаем существующие точки
    const pointsContainer = document.getElementById('points-container');
    pointsContainer.innerHTML = '';
    pointCounter = 0;

    // Добавляем точки из наблюдений
    observations.forEach((obs, index) => {
        pointCounter++;
        const newPoint = document.createElement('div');
        newPoint.className = 'point-row';
        newPoint.innerHTML = `
            <div class="point-label">Точка ${pointCounter}:</div>
            <input type="datetime-local" id="time${pointCounter}" value="${obs.time.replace(' ', 'T')}">
            <input type="number" id="ra${pointCounter}" placeholder="Прямое восхождение (часы)" step="0.1" value="${obs.ra}">
            <input type="number" id="dec${pointCounter}" placeholder="Склонение (градусы)" step="0.1" value="${obs.dec}">
        `;
        pointsContainer.appendChild(newPoint);
    });
}

function fillOrbitResults(orbit) {
    document.getElementById('semiMajorAxis').textContent = orbit.semi_major_axis;
    document.getElementById('eccentricity').textContent = orbit.eccentricity;
    document.getElementById('inclination').textContent = orbit.inclination;
    document.getElementById('longitudeNode').textContent = orbit.longitude_ascending;
    document.getElementById('argumentPerihelion').textContent = orbit.argument_pericenter;
    document.getElementById('trueAnomaly').textContent = orbit.true_anomaly;

    currentOrbitData = orbit;
}

async function savePlanet() {
  const name = document.getElementById('planetName').value.trim();
  if (!name) {
      alert('Введите название планеты');
      return;
  }

  const observations = collectObservationData();
  if (observations.length < 5) {
      alert('Нужно минимум 5 наблюдений для сохранения');
      return;
  }

  if (!currentOrbitData) {
      alert('Сначала рассчитайте параметры орбиты');
      return;
  }

  // ПОЛУЧАЕМ ИЗОБРАЖЕНИЕ ИЗ LOCALSTORAGE
  const imageData = localStorage.getItem('cometImage') || '';

  try {
      const response = await fetch('http://127.0.0.1:5001/api/planets', {
          method: 'POST',
          headers: {
              'Content-Type': 'application/json',
          },
          body: JSON.stringify({
              name: name,
              observations: observations,
              orbital_elements: currentOrbitData,
              image_data: imageData  // ДОБАВЛЯЕМ ИЗОБРАЖЕНИЕ
          })
      });

      const result = await response.json();

      if (result.success) {
          alert('✅ Планета сохранена! ID: ' + result.planet_id);
          document.getElementById('planetName').value = '';
          loadPlanets();
      } else {
          alert('Ошибка сохранения: ' + result.error);
      }
  } catch (error) {
      alert('Ошибка соединения: ' + error.message);
  }
}

// СУЩЕСТВУЮЩИЕ ФУНКЦИИ (оставляем без изменений)

function addPoint() {
    pointCounter++;

    const pointsContainer = document.getElementById("points-container");

    const newPoint = document.createElement("div");
    newPoint.className = "point-row";
    newPoint.innerHTML = `
        <div class="point-label">Точка ${pointCounter}:</div>
        <input type="datetime-local" id="time${pointCounter}">
        <input type="number" id="ra${pointCounter}" placeholder="Прямое восхождение (часы)" step="0.1">
        <input type="number" id="dec${pointCounter}" placeholder="Склонение (градусы)" step="0.1">
    `;

    pointsContainer.appendChild(newPoint);
}

function collectObservationData() {
    const observations = [];

    for (let i = 1; i <= pointCounter; i++) {
        const timeInput = document.getElementById("time" + i);
        const raInput = document.getElementById("ra" + i);
        const decInput = document.getElementById("dec" + i);

        if (!timeInput || !raInput || !decInput) {
            console.warn(`Элементы для точки ${i} не найдены`);
            continue;
        }

        const time = timeInput.value;
        const ra = raInput.value;
        const dec = decInput.value;

        if (time && time.trim() !== '' &&
            ra && ra.trim() !== '' &&
            dec && dec.trim() !== '') {

            observations.push({
                time: time.replace('T', ' ') + ':00',
                ra: parseFloat(ra),
                dec: parseFloat(dec)
            });
        }
    }

    console.log("Собрано наблюдений:", observations.length, observations);
    return observations;
}

async function calculateOrbit() {
    const observations = collectObservationData();
    console.log("Отправляемые данные:", observations);

    if (observations.length < 5) {
        alert('Нужно минимум 5 наблюдений! Заполнено: ' + observations.length);
        return;
    }

    try {
        const response = await fetch('http://127.0.0.1:5001/api/calculate-orbit', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                observations: observations
            })
        });

        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }

        const result = await response.json();
        console.log("Ответ от сервера:", result);

        if (result.success) {
            // Обновляем основные параметры орбиты
            document.getElementById('semiMajorAxis').textContent = result.orbit.semi_major_axis?.toFixed(6) || '-';
            document.getElementById('eccentricity').textContent = result.orbit.eccentricity?.toFixed(6) || '-';
            document.getElementById('inclination').textContent = result.orbit.inclination?.toFixed(6) || '-';
            document.getElementById('longitudeNode').textContent = result.orbit.longitude_ascending?.toFixed(6) || '-';
            document.getElementById('argumentPerihelion').textContent = result.orbit.argument_pericenter?.toFixed(6) || '-';
            document.getElementById('trueAnomaly').textContent = result.orbit.true_anomaly?.toFixed(6) || '-';

            currentOrbitData = result.orbit;
            alert('✅ Орбитальные параметры успешно рассчитаны!');
        } else {
            alert('Ошибка сервера: ' + result.error);
        }
    } catch (error) {
        console.error("Полная ошибка:", error);
        alert('Ошибка соединения: ' + error.message);
    }
}

async function calculateApproach() {
    const semiMajorAxis = document.getElementById('semiMajorAxis').textContent;
    const eccentricity = document.getElementById('eccentricity').textContent;

    if (semiMajorAxis === '-' || eccentricity === '-') {
        alert('Сначала рассчитайте параметры орбиты!');
        return;
    }

    const orbitParams = {
        semi_major_axis: parseFloat(semiMajorAxis),
        eccentricity: parseFloat(eccentricity),
        inclination: parseFloat(document.getElementById('inclination').textContent),
        longitude_ascending: parseFloat(document.getElementById('longitudeNode').textContent),
        argument_pericenter: parseFloat(document.getElementById('argumentPerihelion').textContent),
        true_anomaly: parseFloat(document.getElementById('trueAnomaly').textContent) || 0
    };

    try {
        const response = await fetch('http://127.0.0.1:5001/api/calculate-approach', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                orbit: orbitParams
            })
        });

        const result = await response.json();

        if (result.success) {
            document.getElementById('approachDate').textContent = result.approach.date;
            document.getElementById('approachDistance').textContent = result.approach.distance_au?.toFixed(6) + ' а.е.';
            document.getElementById('collisionStatus').textContent = result.approach.is_safe ? 'Безопасно' : 'Опасно!';
            document.getElementById('collisionStatus').className = result.approach.is_safe ? 'safe-status' : 'danger-status';

            alert('✅ Сближение с Землей рассчитано!');
        } else {
            alert('Ошибка: ' + result.error);
        }
    } catch (error) {
        alert('Ошибка соединения с сервером: ' + error.message);
    }
}



document.addEventListener('DOMContentLoaded', function() {
  // Создаем 5 пустых точек при загрузке
  for (let i = 4; i <= 5; i++) {
      addPoint();
  }

  // Инициализируем загрузку изображений
  initImageUpload();

  // Восстанавливаем сохраненное изображение
  restoreImage();

  loadPlanets(); // Загружаем список планет при запуске
});
