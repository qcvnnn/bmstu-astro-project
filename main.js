let pointCounter = 5;
let currentOrbitData = null;
let currentCometId = null;
let cometsListVisible = false;

// ФУНКЦИИ ДЛЯ УПРАВЛЕНИЯ КОМЕТАМИ

function toggleCometsList() {
    const cometsList = document.getElementById('cometsList');
    cometsListVisible = !cometsListVisible;

    if (cometsListVisible) {
        cometsList.style.display = 'block';
        loadComets();
    } else {
        cometsList.style.display = 'none';
    }
}

async function loadComets() {
    try {
        const response = await fetch('http://127.0.0.1:5001/api/planets');
        const result = await response.json();

        if (result.success) {
            displayComets(result.planets);
        } else {
            showNotification('❌ Ошибка загрузки комет: ' + result.error, 'error');
        }
    } catch (error) {
        showNotification('❌ Ошибка соединения: ' + error.message, 'error');
    }
}

function displayComets(comets) {
    const cometsContainer = document.getElementById('comets-list-container');
    cometsContainer.innerHTML = '';

    if (comets.length === 0) {
        cometsContainer.innerHTML = '<p style="text-align: center; color: #666; padding: 20px;">Нет сохраненных комет</p>';
        return;
    }

    comets.forEach(comet => {
        const cometElement = document.createElement('div');
        cometElement.className = 'comet-card';

        const imageHtml = comet.image_data ?
            `<img src="${comet.image_data}" alt="${comet.name}" class="comet-image">` :
            '<div class="comet-no-image">🌠</div>';

        cometElement.innerHTML = `
            ${imageHtml}
            <div class="comet-info">
                <div class="comet-name">${comet.name}</div>
                <div class="comet-details">
                    Наблюдений: ${comet.observations.length} |
                    Создана: ${new Date(comet.created_at).toLocaleDateString()}
                </div>
            </div>
            <div class="comet-actions">
                <button class="load-comet-btn" onclick="loadCometData(${comet.id})">📊 Загрузить</button>
                <button class="delete-comet-btn" onclick="deleteComet(${comet.id})">🗑️ Удалить</button>
            </div>
        `;

        cometsContainer.appendChild(cometElement);
    });
}

async function loadCometData(cometId) {
    try {
        const response = await fetch('http://127.0.0.1:5001/api/planets');
        const result = await response.json();

        if (result.success) {
            const comet = result.planets.find(p => p.id === cometId);
            if (comet) {
                // Заполняем данные кометы в поля формы
                fillObservations(comet.observations);
                fillOrbitResults(comet.orbital_elements);

                // Загружаем изображение
                if (comet.image_data) {
                    loadCometImage(comet.image_data);
                } else {
                    removeImage();
                }

                currentCometId = cometId;

                showNotification(`✅ Данные кометы "${comet.name}" загружены`, 'success');
                // НЕ закрываем список комет
            }
        }
    } catch (error) {
        showNotification('❌ Ошибка загрузки данных: ' + error.message, 'error');
    }
}

async function deleteComet(cometId) {
    if (!confirm('Удалить эту комету?')) return;

    try {
        const response = await fetch(`http://127.0.0.1:5001/api/planets/${cometId}`, {
            method: 'DELETE'
        });
        const result = await response.json();

        if (result.success) {
            showNotification('🗑️ Комета удалена', 'success');
            // Просто обновляем список, не закрывая его
            loadComets();
        } else {
            showNotification('❌ Ошибка удаления: ' + result.error, 'error');
        }
    } catch (error) {
        showNotification('❌ Ошибка соединения: ' + error.message, 'error');
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

function loadCometImage(imageData) {
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

async function saveCometFinal() {
    const cometName = prompt('Введите название кометы:');
    if (!cometName || !cometName.trim()) {
        showNotification('❌ Название кометы обязательно', 'error');
        return;
    }

    // Проверяем, существует ли уже комета с таким названием
    const isDuplicate = await checkDuplicateComet(cometName.trim());
    if (isDuplicate) {
        showNotification('❌ Комета с таким названием уже существует', 'error');
        return;
    }

    await saveCometToDatabase(cometName.trim());
}

async function checkDuplicateComet(cometName) {
    try {
        const response = await fetch('http://127.0.0.1:5001/api/planets');
        const result = await response.json();

        if (result.success) {
            const existingComet = result.planets.find(comet =>
                comet.name.toLowerCase() === cometName.toLowerCase()
            );
            return !!existingComet;
        }
        return false;
    } catch (error) {
        console.error('Ошибка проверки дубликата:', error);
        return false;
    }
}

async function saveCometToDatabase(cometName) {
    const observations = collectObservationData();
    if (observations.length < 5) {
        showNotification('❌ Нужно минимум 5 наблюдений для сохранения', 'error');
        return;
    }

    if (!currentOrbitData) {
        showNotification('❌ Сначала рассчитайте параметры орбиты', 'error');
        return;
    }

    // Получаем изображение из localStorage
    const imageData = localStorage.getItem('cometImage') || '';

    try {
        const response = await fetch('http://127.0.0.1:5001/api/planets', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                name: cometName,
                observations: observations,
                orbital_elements: currentOrbitData,
                image_data: imageData
            })
        });

        const result = await response.json();

        if (result.success) {
            showNotification(`✅ Комета "${cometName}" сохранена! ID: ${result.planet_id}`, 'success');
            loadComets(); // Обновляем список комет
        } else {
            showNotification('❌ Ошибка сохранения: ' + result.error, 'error');
        }
    } catch (error) {
        showNotification('❌ Ошибка соединения: ' + error.message, 'error');
    }
}

// ФУНКЦИИ ДЛЯ РАБОТЫ С ИЗОБРАЖЕНИЯМИ
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
            <div class="placeholder-icon">🛸</div>
            <p>Перетащите сюда изображение<br>или нажмите кнопку ниже</p>
        </div>
    `;
    imageUpload.value = '';

    // Удаляем из localStorage
    localStorage.removeItem('cometImage');

    showNotification('🗑️ Изображение удалено', 'info');
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

// ОСНОВНЫЕ ФУНКЦИИ ПРИЛОЖЕНИЯ
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
            showNotification('✅ Орбитальные параметры успешно рассчитаны!', 'success');
        } else {
            showNotification('❌ Ошибка сервера: ' + result.error, 'error');
        }
    } catch (error) {
        console.error("Полная ошибка:", error);
        showNotification('❌ Ошибка соединения: ' + error.message, 'error');
    }
}

async function calculateApproach() {
    const semiMajorAxis = document.getElementById('semiMajorAxis').textContent;
    const eccentricity = document.getElementById('eccentricity').textContent;

    if (semiMajorAxis === '-' || eccentricity === '-') {
        showNotification('❌ Сначала рассчитайте параметры орбиты!', 'error');
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

            showNotification('✅ Сближение с Землей рассчитано!', 'success');
        } else {
            showNotification('❌ Ошибка: ' + result.error, 'error');
        }
    } catch (error) {
        showNotification('❌ Ошибка соединения с сервером: ' + error.message, 'error');
    }
}

// ИНИЦИАЛИЗАЦИЯ ПРИ ЗАГРУЗКЕ СТРАНИЦЫ
document.addEventListener('DOMContentLoaded', function() {
    // Инициализируем загрузку изображений
    initImageUpload();

    // Восстанавливаем сохраненное изображение
    restoreImage();

    // Загружаем список комет при запуске (но не показываем)
    loadComets();
});
