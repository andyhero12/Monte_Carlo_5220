import subprocess
import numpy as np
#import matplotlib.pyplot as plt

# Define the range of particle sizes to test
particle_sizes = [1e3, 1e6, 2e6]  # Geometrically spaced sizes for log-log plot

# Initialize an empty list to store simulation times
simulation_times = []

# Loop through each particle size, run the simulation, and capture the output time
for size in particle_sizes:
    print(f"Size is {size}")
    command = f"./build/serial -s 1 -n {size}"
    # Adjusted for compatibility with Python < 3.7
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    output = result.stdout
    # print(output)
    # Extract the simulation time from the command output
    time_str = output.split("Simulation Time = ")[1].split(" seconds")[0]
    simulation_time = float(time_str)
    simulation_times.append(simulation_time)
    print(f"Simulation time for {size} particles: {simulation_time} seconds")


# Calculate the slope of the line of best fit on the log-log scale
log_particle_sizes = np.log(particle_sizes)
log_simulation_times = np.log(simulation_times)
slope, intercept = np.polyfit(log_particle_sizes, log_simulation_times, 1)

print(f"Slope of the line of best fit (log-log scale): {slope}")
print(f"Intercept: {intercept}")

# Plotting (commented out)
# plt.figure(figsize=(10, 6))
# plt.loglog(particle_sizes, simulation_times, marker='o', linestyle='-')
# plt.title('Simulation Time vs Number of Particles')
# plt.xlabel('Number of Particles')
# plt.ylabel('Simulation Time (seconds)')
# plt.grid(True, which="both", ls="--")
# plt.savefig("performance_plot.png")
# plt.show()