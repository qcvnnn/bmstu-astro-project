// let pointCounter = 5;
// let currentOrbitData = null;

// async function loadPlanets() {
//   try {
//       const response = await fetch('http://127.0.0.1:5001/api/planets');
//       const result = await response.json();

//       if (result.success) {
//           displayPlanets(result.planets);
//       } else {
//           alert('Ошибка загрузки планет: ' + result.error);
//       }
//   } catch (error) {
//       alert('Ошибка соединения: ' + error.message);
//   }
// }

// function displayPlanets(planets) {
//   const planetsList = document.getElementById('planets-list');
//   planetsList.innerHTML = '';

//   if (planets.length === 0) {
//       planetsList.innerHTML = '<p>Нет сохраненных планет</p>';
//       return;
//   }

//   planets.forEach(planet => {
//       const planetElement = document.createElement('div');
//       planetElement.className = 'planet-card';
//       planetElement.innerHTML = `
//           <div class="planet-header">
//               <h3>${planet.name}</h3>
//               <button onclick="deletePlanet(${planet.id})">🗑️</button>
//           </div>
//           <div class="planet-info">
//               <p><strong>Наблюдения:</strong> ${planet.observations.length} точек</p>
//               <p><strong>Большая полуось:</strong> ${planet.orbital_elements.semi_major_axis} а.е.</p>
//               <p><strong>Эксцентриситет:</strong> ${planet.orbital_elements.eccentricity}</p>
//           </div>
//           <button onclick="loadPlanetData(${planet.id})">📊 Загрузить данные</button>
//       `;
//       planetsList.appendChild(planetElement);
//   });
// }

// async function deletePlanet(planetId) {
//   if (!confirm('Удалить эту планету?')) return;

//   try {
//       const response = await fetch(`http://127.0.0.1:5001/api/planets/${planetId}`, {
//           method: 'DELETE'
//       });
//       const result = await response.json();

//       if (result.success) {
//           alert('Планета удалена');
//           loadPlanets();
//       } else {
//           alert('Ошибка удаления: ' + result.error);
//       }
//   } catch (error) {
//       alert('Ошибка соединения: ' + error.message);
//   }
// }

// async function loadPlanetData(planetId) {
//   try {
//       const response = await fetch('http://127.0.0.1:5001/api/planets');
//       const result = await response.json();

//       if (result.success) {
//           const planet = result.planets.find(p => p.id === planetId);
//           if (planet) {
//               // Заполняем поля наблюдениями
//               fillObservations(planet.observations);
//               // Заполняем результаты орбиты
//               fillOrbitResults(planet.orbital_elements);
//               alert(`Данные планеты "${planet.name}" загружены`);
//           }
//       }
//   } catch (error) {
//       alert('Ошибка загрузки данных: ' + error.message);
//   }
// }

// function fillObservations(observations) {
//   // Очищаем существующие точки
//   const pointsContainer = document.getElementById('points-container');
//   pointsContainer.innerHTML = '';
//   pointCounter = 0;

//   // Добавляем точки из наблюдений
//   observations.forEach((obs, index) => {
//       pointCounter++;
//       const newPoint = document.createElement('div');
//       newPoint.className = 'point-row';
//       newPoint.innerHTML = `
//           <div class="point-label">Точка ${pointCounter}:</div>
//           <input type="datetime-local" id="time${pointCounter}" value="${obs.time.replace(' ', 'T')}">
//           <input type="number" id="ra${pointCounter}" placeholder="Прямое восхождение (часы)" step="0.1" value="${obs.ra}">
//           <input type="number" id="dec${pointCounter}" placeholder="Склонение (градусы)" step="0.1" value="${obs.dec}">
//       `;
//       pointsContainer.appendChild(newPoint);
//   });
// }

// function fillOrbitResults(orbit) {
//   document.getElementById('semiMajorAxis').textContent = orbit.semi_major_axis;
//   document.getElementById('eccentricity').textContent = orbit.eccentricity;
//   document.getElementById('inclination').textContent = orbit.inclination;
//   document.getElementById('longitudeNode').textContent = orbit.longitude_ascending;
//   document.getElementById('argumentPerihelion').textContent = orbit.argument_pericenter;

//   currentOrbitData = orbit;
// }

// async function savePlanet() {
//   const name = document.getElementById('planetName').value.trim();
//   if (!name) {
//       alert('Введите название планеты');
//       return;
//   }

//   const observations = collectObservationData();
//   if (observations.length < 3) {
//       alert('Нужно минимум 3 наблюдения для сохранения');
//       return;
//   }

//   if (!currentOrbitData) {
//       alert('Сначала рассчитайте параметры орбиты');
//       return;
//   }

//   try {
//       const response = await fetch('http://127.0.0.1:5001/api/planets', {
//           method: 'POST',
//           headers: {
//               'Content-Type': 'application/json',
//           },
//           body: JSON.stringify({
//               name: name,
//               observations: observations,
//               orbital_elements: currentOrbitData
//           })
//       });

//       const result = await response.json();

//       if (result.success) {
//           alert('Планета сохранена! ID: ' + result.planet_id);
//           document.getElementById('planetName').value = '';
//           loadPlanets();
//       } else {
//           alert('Ошибка сохранения: ' + result.error);
//       }
//   } catch (error) {
//       alert('Ошибка соединения: ' + error.message);
//   }
// }

// function addPoint() {
//   pointCounter++;

//   const pointsContainer = document.getElementById("points-container");

//   const newPoint = document.createElement("div");
//   newPoint.className = "point-row";
//   newPoint.innerHTML = `
//         <div class="point-label">Точка ${pointCounter}:</div>
//         <input type="datetime-local" id="time${pointCounter}">
//         <input type="number" id="ra${pointCounter}" placeholder="Прямое восхождение (часы)" step="0.1">
//         <input type="number" id="dec${pointCounter}" placeholder="Склонение (градусы)" step="0.1">
//     `;

//   pointsContainer.appendChild(newPoint);
// }

// function collectObservationData() {
//   const observations = [];

//   for (let i = 1; i <= pointCounter; i++) {
//       const timeInput = document.getElementById("time" + i);
//       const raInput = document.getElementById("ra" + i);
//       const decInput = document.getElementById("dec" + i);

//       // ПРОВЕРЯЕМ СУЩЕСТВОВАНИЕ ЭЛЕМЕНТОВ
//       if (!timeInput || !raInput || !decInput) {
//           console.warn(`Элементы для точки ${i} не найдены`);
//           continue;
//       }

//       const time = timeInput.value;
//       const ra = raInput.value;
//       const dec = decInput.value;

//       // БОЛЕЕ СТРОГАЯ ПРОВЕРКА
//       if (time && time.trim() !== '' &&
//           ra && ra.trim() !== '' &&
//           dec && dec.trim() !== '') {

//           observations.push({
//               time: time.replace('T', ' ') + ':00',
//               ra: parseFloat(ra),
//               dec: parseFloat(dec)
//           });
//       }
//   }

//   console.log("Собрано наблюдений:", observations.length, observations);
//   return observations;
// }

// async function calculateOrbit() {
//   const observations = collectObservationData();
//   console.log("Отправляемые данные:", observations); // ДЛЯ ОТЛАДКИ

//   if (observations.length < 3) {
//       alert('Нужно минимум 3 наблюдения! Заполнено: ' + observations.length);
//       return;
//   }

//   try {
//       const response = await fetch('http://127.0.0.1:5001/api/calculate-orbit', {
//           method: 'POST',
//           headers: {
//               'Content-Type': 'application/json',
//           },
//           body: JSON.stringify({
//               observations: observations
//           })
//       });

//       // ПРОВЕРЯЕМ STATUS RESPONSE
//       if (!response.ok) {
//           throw new Error(`HTTP error! status: ${response.status}`);
//       }

//       const result = await response.json();
//       console.log("Ответ от сервера:", result);

//       if (result.success) {
//           document.getElementById('semiMajorAxis').textContent = result.orbit.semi_major_axis;
//           document.getElementById('eccentricity').textContent = result.orbit.eccentricity;
//           document.getElementById('inclination').textContent = result.orbit.inclination;
//           document.getElementById('longitudeNode').textContent = result.orbit.longitude_ascending;
//           document.getElementById('argumentPerihelion').textContent = result.orbit.argument_pericenter;

//           alert('✅ Орбитальные параметры успешно рассчитаны!');
//       } else {
//           alert('Ошибка сервера: ' + result.error);
//       }
//   } catch (error) {
//       console.error("Полная ошибка:", error);
//       alert('Ошибка соединения: ' + error.message);
//   }
// }

// async function calculateApproach() {
//     const semiMajorAxis = document.getElementById('semiMajorAxis').textContent;
//     const eccentricity = document.getElementById('eccentricity').textContent;
//     const inclination = document.getElementById('inclination').textContent;
//     const longitudeNode = document.getElementById('longitudeNode').textContent;
//     const argumentPerihelion = document.getElementById('argumentPerihelion').textContent;

//     if (semiMajorAxis === '-' || eccentricity === '-') {
//         alert('Сначала рассчитайте параметры орбиты!');
//         return;
//     }

//     const orbitParams = {
//         semi_major_axis: parseFloat(semiMajorAxis),
//         eccentricity: parseFloat(eccentricity),
//         inclination: parseFloat(inclination),
//         longitude_ascending: parseFloat(longitudeNode),
//         argument_pericenter: parseFloat(argumentPerihelion)
//     };

//     try {
//         const response = await fetch('http://127.0.0.1:5001/api/calculate-approach', {
//             method: 'POST',
//             headers: {
//                 'Content-Type': 'application/json',
//             },
//             body: JSON.stringify({
//                 orbit: orbitParams
//             })
//         });

//         const result = await response.json();

//         if (result.success) {
//             document.getElementById('approachDate').textContent = result.approach.date;
//             document.getElementById('approachDistance').textContent = result.approach.distance_au + ' а.е.';
//             document.getElementById('collisionStatus').textContent = result.approach.is_safe ? 'Безопасно' : 'Опасно!';
//             document.getElementById('collisionStatus').className = result.approach.is_safe ? 'safe-status' : 'danger-status';

//             alert('✅ Сближение с Землей рассчитано!');
//         } else {
//             alert('Ошибка: ' + result.error);
//         }
//     } catch (error) {
//         alert('Ошибка соединения с сервером: ' + error.message);
//     }
// }

// // Старая функция (можно удалить или оставить как альтернативу)
// function sendToServer() {
//   const observations = collectObservationData();

//   if (observations.length < 3) {
//     alert('Нужно минимум 3 наблюдения!');
//     return;
//   }

//   fetch('http://127.0.0.1:5001/api/calculate-orbit', {
//     method: 'POST',
//     headers: {
//       'Content-Type': 'application/json',
//     },
//     body: JSON.stringify({ observations: observations })
//   })
//   .then(response => response.json())
//   .then(result => {
//     if (result.success) {
//       document.getElementById('semiMajorAxis').textContent = result.orbit.semi_major_axis;
//       document.getElementById('eccentricity').textContent = result.orbit.eccentricity;
//       document.getElementById('inclination').textContent = result.orbit.inclination;
//       document.getElementById('longitudeNode').textContent = result.orbit.longitude_ascending;
//       document.getElementById('argumentPerihelion').textContent = result.orbit.argument_pericenter;
//     } else {
//       alert('Ошибка: ' + result.error);
//     }
//   })
//   .catch(error => {
//     alert('Ошибка соединения с сервером: ' + error.message);
//   });
// }

// document.addEventListener('DOMContentLoaded', function() {
//     setTimeout(() => {
//         fillTestData();
//     }, 1000);
// });

// function fillTestData() {
//     const testData = [
//         { time: '2025-10-25T00:00', ra: '15.30977', dec: '-18.61633' },
//         { time: '2025-10-27T00:00', ra: '15.40572', dec: '-18.99403' },
//         { time: '2025-10-29T00:00', ra: '15.50238', dec: '-19.36158' },
//         { time: '2025-10-31T00:00', ra: '15.59917', dec: '-19.71861' },
//         { time: '2025-11-02T00:00', ra: '15.86444', dec: '-20.06417' }
//     ];

//     for (let i = 0; i < testData.length; i++) {
//         const point = testData[i];
//         if (document.getElementById('time' + (i + 1))) {
//             document.getElementById('time' + (i + 1)).value = point.time;
//             document.getElementById('ra' + (i + 1)).value = point.ra;
//             document.getElementById('dec' + (i + 1)).value = point.dec;
//         }
//     }

//     console.log('✅ Тестовые данные загружены');
// }

let pointCounter = 5;

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

    if (observations.length < 3) {
        alert('Нужно минимум 3 наблюдения! Заполнено: ' + observations.length);
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
    const inclination = document.getElementById('inclination').textContent;
    const longitudeNode = document.getElementById('longitudeNode').textContent;
    const argumentPerihelion = document.getElementById('argumentPerihelion').textContent;
    const trueAnomaly = document.getElementById('trueAnomaly').textContent;    if (semiMajorAxis === '-' || eccentricity === '-') {
      alert('Сначала рассчитайте параметры орбиты!');
      return;
  }

  const orbitParams = {
      semi_major_axis: parseFloat(semiMajorAxis),
      eccentricity: parseFloat(eccentricity),
      inclination: parseFloat(inclination),
      longitude_ascending: parseFloat(longitudeNode),
      argument_pericenter: parseFloat(argumentPerihelion),
      true_anomaly: parseFloat(trueAnomaly) || 0
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

// Старая функция для совместимости
function sendToServer() {
  calculateOrbit();
}

// АВТОМАТИЧЕСКОЕ ЗАПОЛНЕНИЕ ТЕСТОВЫМИ ДАННЫМИ
document.addEventListener('DOMContentLoaded', function() {
  setTimeout(() => {
      fillTestData();
  }, 1000);
});
