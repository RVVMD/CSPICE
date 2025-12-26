import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation, PillowWriter

# ПАРАМЕТРЫ
data_file = "results.csv"
output_gif = "wave_propagation.gif"
fps = 30
line_length_km = 100  # длина линии в км
max_voltage_kV = 1.75  # максимальное напряжение для масштаба оси Y
every_nth_frame = 1   # использовать каждый N-й кадр для ускорения

# Загрузка данных
df = pd.read_csv(data_file)
times = df["Time(us)"].values[::every_nth_frame]
voltages = df.filter(like="V_node_").values[::every_nth_frame] / 1000  # в кВ

# Настройка графика
fig, ax = plt.subplots(figsize=(12, 7))
line, = ax.plot([], [], 'b-', linewidth=2.5)
ax.set_xlim(0, line_length_km)
ax.set_ylim(-max_voltage_kV, max_voltage_kV)
ax.set_xlabel('Расстояние по линии (км)', fontsize=12)
ax.set_ylabel('Напряжение (кВ)', fontsize=12)
ax.set_title('Распространение волны при ударе молнии', fontsize=14)
ax.grid(True, linestyle='--', alpha=0.7)
ax.axvline(x=line_length_km/2.5, color='r', linestyle='--', alpha=0.3, label='Точка удара')
ax.legend()

# Инициализация
def init():
    line.set_data([], [])
    return line,

# Обновление кадра
def update(frame):
    x = np.linspace(0, line_length_km, voltages.shape[1])
    y = voltages[frame]

    line.set_data(x, y)
    ax.set_title(f'Распространение волны при ударе молнии\nВремя: {times[frame]:.2f} мкс', fontsize=14)
    return line,

# Создание анимации
ani = FuncAnimation(fig, update, frames=len(times), init_func=init,
                    blit=True, interval=1000/fps)

# Сохранение GIF
writer = PillowWriter(fps=fps)
ani.save(output_gif, writer=writer, dpi=100)

print(f"Анимация сохранена в {output_gif}")
print(f"Количество кадров: {len(times)}")
