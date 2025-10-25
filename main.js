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

  function collectObservationData() {
    const observations = [];

    // Собираем данные со всех точек (от 1 до pointCounter)
    for (let i = 1; i <= pointCounter; i++) {
      const time = document.getElementById("time" + i)?.value;
      const ra = document.getElementById("ra" + i)?.value;
      const dec = document.getElementById("dec" + i)?.value;

      // Добавляем только если есть данные
      if (time && ra && dec) {
        observations.push({
          time: time,
          ra: ra,
          dec: dec,
        });
      }
    }
    return observations;
  }
}
