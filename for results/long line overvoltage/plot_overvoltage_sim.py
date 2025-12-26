import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import imageio.v2 as imageio

# Настройка стиля matplotlib
mpl.rcParams['font.family'] = 'DejaVu Sans'
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14

# Загрузка данных анимации
data = pd.read_csv('animation_data.csv')

# Извлечение параметров
num_frames = int(data['frame'].max()) + 1
num_segments = (data.shape[1] - 2)
line_length = 3000  # км

# Создание массива расстояний для каждого сегмента
distances = np.linspace(0, line_length, num_segments)

# Создание фигуры и осей
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), height_ratios=[3, 1], gridspec_kw={'hspace': 0.3})
fig.suptitle('Перенапряжение на конце линии при отключении нагрузки', fontsize=14)

# Основной график напряжения вдоль линии
line, = ax1.plot([], [], 'b-', linewidth=2.0)
ax1.set_xlim(0, line_length)
ax1.set_ylim(-600000, 600000)
ax1.set_xlabel('Расстояние от начала линии (км)', fontsize=12)
ax1.set_ylabel('Напряжение (В)', fontsize=12)
ax1.grid(True, linestyle='--', alpha=0.7)

# График напряжения в конце линии во времени
time_data = []
end_voltage_data = []
end_line, = ax2.plot([], [], 'g-', linewidth=2.0)
ax2.set_xlim(0, 0.15)
ax2.set_ylim(-600000, 600000)
ax2.set_xlabel('Время (с)', fontsize=12)
ax2.set_ylabel('Напряжение в конце (В)', fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.7)

# Вертикальная линия на 0.1с (момент отключения нагрузки)
ax2.axvline(x=0.1, color='r', linestyle='--', alpha=0.7, label='Отключение нагрузки')
ax2.legend(loc='upper right')

# Индикатор максимального напряжения
max_voltage_text = ax2.text(0.02, 0.85, '', transform=ax2.transAxes, fontsize=10,
                           bbox=dict(facecolor='white', alpha=0.8))
max_voltage_line = ax2.axhline(y=0, color='purple', linestyle=':', alpha=0.7)

# Текстовые аннотации
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes, fontsize=10,
                    bbox=dict(facecolor='white', alpha=0.8))

# Список для хранения кадров GIF
frames = []
max_end_voltage = 0
time_points = data['time'].unique()
time_points.sort()

# Создание всех кадров
for frame_idx in range(num_frames):
    frame_data = data[data['frame'] == frame_idx]
    if frame_data.empty or len(frame_data) == 0:
        continue

    time_val = frame_data['time'].values[0]
    voltages = frame_data.iloc[:, 2:].values[0]

    # Обновление данных линии
    line.set_data(distances, voltages)

    # Обновление времени
    time_text.set_text(f'Время: {time_val*1000:.2f} мс')

    # Обновление данных для графика напряжения в конце линии
    end_voltage = voltages[-1]
    time_data.append(time_val)
    end_voltage_data.append(end_voltage)

    # Обновление максимального напряжения в конце линии
    if abs(end_voltage) > abs(max_end_voltage):
        max_end_voltage = end_voltage

    end_line.set_data(time_data, end_voltage_data)
    max_voltage_text.set_text(f'Макс. напряжение в конце: {abs(max_end_voltage)/1000:.1f} кВ')
    max_voltage_line.set_ydata([max_end_voltage])

    # Рендеринг кадра
    fig.canvas.draw()

    # Преобразование в изображение
    image = np.array(fig.canvas.buffer_rgba())
    image = image[:, :, :3]
    frames.append(image.copy())

    # Обновление прогресса каждые 1%
    if frame_idx % max(1, num_frames // 100) == 0:
        percent = int(frame_idx * 100 / num_frames)
        print(f"\r{percent}%", end='', flush=True)

print("\r100%")

# Сохранение GIF файла
imageio.mimsave('transmission_line_animation.gif', frames, fps=25, loop=0)
