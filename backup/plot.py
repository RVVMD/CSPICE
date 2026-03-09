import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import os

# Configuration based on TR5_sim.c
LINE1_SEGMENTS = 400
LINE2_SEGMENTS = 240
LINE3_SEGMENTS = 120

LINE1_LENGTH_KM = 100.0
LINE2_LENGTH_KM = 60.0
LINE3_LENGTH_KM = 30.0

CSV_FILE = "voltage_results.csv"

def main():
    if not os.path.exists(CSV_FILE):
        print(f"Error: {CSV_FILE} not found. Please run the C simulation first.")
        return

    # Load data
    df = pd.read_csv(CSV_FILE)
    data = df.to_numpy()

    # Time column is index 0
    time_data = data[:, 0]

    # Define spatial coordinates for each line
    # Line 1: 0 to 100 km
    x1 = np.linspace(0, LINE1_LENGTH_KM, LINE1_SEGMENTS + 1)
    # Line 2: 100 to 160 km
    x2 = np.linspace(LINE1_LENGTH_KM, LINE1_LENGTH_KM + LINE2_LENGTH_KM, LINE2_SEGMENTS + 1)
    # Line 3: 160 to 190 km
    x3 = np.linspace(LINE1_LENGTH_KM + LINE2_LENGTH_KM,
                     LINE1_LENGTH_KM + LINE2_LENGTH_KM + LINE3_LENGTH_KM,
                     LINE3_SEGMENTS + 1)

    # Combine spatial coordinates for continuous plotting
    x_total = np.concatenate([x1, x2, x3])

    # Column indices for voltage data
    # Line 1: 1 to 401 (401 nodes)
    idx1_start = 1
    idx1_end = 1 + LINE1_SEGMENTS + 1
    # Line 2: 402 to 642 (241 nodes)
    idx2_start = idx1_end
    idx2_end = idx2_start + LINE2_SEGMENTS + 1
    # Line 3: 643 to 763 (121 nodes)
    idx3_start = idx2_end
    idx3_end = idx3_start + LINE3_SEGMENTS + 1

    # Extract voltage matrices
    v1 = data[:, idx1_start:idx1_end]
    v2 = data[:, idx2_start:idx2_end]
    v3 = data[:, idx3_start:idx3_end]

    # Combine voltage data
    v_total = np.concatenate([v1, v2, v3], axis=1)

    # Setup Plot
    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=2)
    ax.set_xlim(0, LINE1_LENGTH_KM + LINE2_LENGTH_KM + LINE3_LENGTH_KM)
    ax.set_ylim(-5000, 40000) # Adjust based on expected voltage range (Source is 35kV)
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Voltage (V)")
    ax.set_title("Transmission Line Voltage Propagation")
    ax.grid(True)

    # Text annotation for time
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

    def init():
        line.set_data([], [])
        time_text.set_text("")
        return line, time_text

    def update(frame):
        # frame corresponds to row index in data
        y = v_total[frame, :]
        line.set_data(x_total, y)
        t = time_data[frame]
        time_text.set_text(f"Time: {t*1000:.3f} ms")
        return line, time_text

    # Create Animation
    # Interval in milliseconds. dt in C is 1e-6, but plotting every step may be too fast.
    # Downsample if necessary. Here we assume every step is plotted.
    anim = animation.FuncAnimation(fig, update, frames=len(time_data),
                                   init_func=init, blit=True, interval=10)

    # Save or Show
    # To save: anim.save('voltage_propagation.mp4', writer='ffmpeg', fps=30)
    plt.show()

if __name__ == "__main__":
    main()
