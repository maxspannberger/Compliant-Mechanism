import numpy as np
from Positional_Analysis import positional_analysis
from Velocity_Analysis import velocity_analysis



t=np.linspace(0, 2 * np.pi, 100)  # Time array for one full rotation
theta2lst = np.linspace(0, 2 * np.pi, 100)  # Input angle array (in radians)
i=0  # Index for the current time step

for t[i] in t:
    theta2 = theta2lst[i]  # Input angle (in radians)
    
    # Define input and output positions based on the current input angle
    pos_in = [np.cos(theta2), np.sin(theta2)]  # Input position (x, y)
    pos_out = [np.cos(theta2 + np.pi / 4), np.sin(theta2 + np.pi / 4)]  # Output position (x, y)
    
    # Perform positional analysis
    theta3, theta4 = positional_analysis(theta2, pos_in, pos_out)
    
    # Perform velocity analysis
    omega2 = 1.0  # Angular velocity of the crank (rad/s)
    omega3, omega4 = velocity_analysis(theta2, theta3, theta4, omega2, pos_out)
    
    # Print results for the current time step
    print(f"Time: {t[i]:.2f} s, Input Angle: {np.degrees(theta2):.2f} degrees, "
          f"Output Angles: theta3 = {np.degrees(theta3):.2f} degrees, "
          f"theta4 = {np.degrees(theta4):.2f} degrees, "
          f"Angular Velocities: omega3 = {omega3:.2f} rad/s, omega4 = {omega4:.2f} rad/s")
    
    i += 1  # Increment index for the next time step