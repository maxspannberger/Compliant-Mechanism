import numpy as np
import math as mt



def positional_analysis(theta2, pos_in, pos_out):
    """
    Function to perform positional analysis of a four-bar linkage mechanism.
    
    Parameters:
    theta2 : float
        Input angle (in radians) of the crank.
    pos_in : list or np.array
        Coordinates of the input point (x, y).
    pos_out : list or np.array
        Coordinates of the output point (x, y).
    
    Returns:
    theta3 : float
        Output angle (in radians) of the coupler.
    theta4 : float
        Output angle (in radians) of the follower.
    """
    
    # Extract coordinates from input and output positions
    x_in, y_in = pos_in
    x_out, y_out = pos_out
    
    # Calculate lengths of the links based on positions
    L1 = np.sqrt((x_in - 0)**2 + (y_in - 0)**2)  # Length of crank
    L2 = np.sqrt((x_out - x_in)**2 + (y_out - y_in)**2)  # Length of coupler
    L3 = np.sqrt((x_out - 0)**2 + (y_out - 0)**2)  # Length of follower
    
    # Calculate angles using the law of cosines
    cos_theta3 = (L1**2 + L2**2 - L3**2) / (2 * L1 * L2)
    theta3 = mt.acos(cos_theta3)
    
    cos_theta4 = (L1**2 + L3**2 - L2**2) / (2 * L1 * L3)
    theta4 = mt.acos(cos_theta4)
    
    return theta3, theta4

if __name__ == "__main__":
    # Example usage
    theta2 = np.radians(30)  # Input angle in radians
    pos_in = [1.0, 0.0]      # Input position (x, y)
    pos_out = [2.0, 1.0]     # Output position (x, y)
    
    theta3, theta4 = positional_analysis(theta2, pos_in, pos_out)
    print(f"Output angles: theta3 = {np.degrees(theta3):.2f} degrees, theta4 = {np.degrees(theta4):.2f} degrees")