import numpy as np

flex_length=0.09 #flexure length in m
m_2=0.05 #mass of the payload in kg, assumed for now
h=0.3 #height of the jump in m
g=9.81 #acceleration due to gravity in m/s^2
init_angle=[10,15,20,25] #assumed inital deflection of the flexure in degrees
density=8000 #density of the flexure material in kg/m^3 - FIND THE VALUES FOR STEEL

thickness = [0.0005, 0.001, 0.0015] #thickness of the flexure in m
width = [0.005, 0.01, 0.015] #width of the flexure in m

dict_flexure = {}

# Add one dictionary entry for every thickness and width combination.
for tick in thickness:
    for wid in width:
        for angle in init_angle:
            m_1=density*tick*wid*flex_length
            stored_energy=h*(m_1+m_2)**2*g/(m_2*np.sin(np.radians(angle))**2)
            dict_flexure[(wid, tick, angle)] = {
                'flex_length': flex_length,
                'init_angle': angle,
                'thickness': tick,
                'width': wid,
                'energy_stored': stored_energy,
            }

print(dict_flexure)