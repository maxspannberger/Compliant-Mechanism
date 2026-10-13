import numpy as np



def velocity_analysis(theta2, theta3, theta4, pos_in, pos_out):
    """
    Function to perform velocity analysis of a four-bar linkage mechanism.
    
    Parameters:
    theta2 : float
        Input angle (in radians) of the crank.
    theta3 : float
        Output angle (in radians) of the coupler.
    theta4 : float
        Output angle (in radians) of the follower.
    pos_in : list or np.array
        Coordinates of the input point (x, y).
    pos_out : list or np.array
        Coordinates of the output point (x, y).
    
    Returns:
    omega3 : float
        Angular velocity (in radians per second) of the coupler.
    omega4 : float
        Angular velocity (in radians per second) of the follower.
    """
    
    # Extract coordinates from input and output positions
    x_in, y_in = pos_in
    x_out, y_out = pos_out
    
    # Calculate lengths of the links based on positions
    L1 = np.sqrt((x_in - 0)**2 + (y_in - 0)**2)  # Length of crank
    L2 = np.sqrt((x_out - x_in)**2 + (y_out - y_in)**2)  # Length of coupler
    L3 = np.sqrt((x_out - 0)**2 + (y_out - 0)**2)  # Length of follower
    
    # Calculate angular velocities using the law of sines
    sin_theta3 = np.sin(theta3)
    sin_theta4 = np.sin(theta4)
    
    omega3 = (L1 * sin_theta3) / (L2 * sin_theta4)
    omega4 = (L1 * sin_theta4) / (L3 * sin_theta3)
    
    return omega3, omega4