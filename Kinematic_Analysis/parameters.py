import numpy as np


def parameters():

    tot_flex_length=0.09 #flexure length in m
    flex_length_initial=tot_flex_length/2 #initial flexure length in m
    m_2=0.05 #mass of the payload in kg, assumed for now
    h=0.3 #height of the jump in m
    g=9.81 #acceleration due to gravity in m/s^2
    gamma=0.85 #fixed-guided boundary condition
    K_Theta=2.65    #fixed-guided boundary condition
    init_angle=[10,15,20,25] #assumed inital deflection of the flexure in degrees
    density=8000 #density of the flexure material in kg/m^3 - FIND THE VALUES FOR STEEL

    thickness = [0.0005, 0.001, 0.0015] #thickness of the flexure in m
    width = [0.005, 0.01, 0.015] #width of the flexure in m

    dict_flexure = {}

    # Add one dictionary entry for every thickness and width combination.
    for tick in thickness:
        for wid in width:
            for angle in init_angle:
                m_1=density*tick*wid*flex_length_initial #mass of the flexure in kg
                stored_energy=h*(m_1+m_2)**2*g/(m_2*np.sin(np.radians(angle))**2)
                flex_length_deformed=flex_length_initial*np.cos(np.radians(angle))
                d=flex_length_deformed-flex_length_initial
                A=wid*tick
                I=wid*tick**3/12
                E=205*10**9 #Pa, Young's modulus of spring steel apparently from 200 to 210 GPa
                K=1 #fixed-guided boundary condition
                P_cr_buckling_preload=np.pi**2*E*I/((K*tot_flex_length)**2)     #critical buckling load of the straight flexure in N, based on the Euler buckling formula
                preload_force=P_cr_buckling_preload*(1+np.radians(angle)**2/8) #post buckling force relation according to Bazant
                dict_flexure[(wid, tick, angle)] = {
                    'flex_length_initial': flex_length_deformed,
                    'init_angle': np.radians(angle),
                    'gamma': gamma,
                    'K_Theta': K_Theta,
                    'thickness': tick,
                    'width': wid,
                    'preload_force': preload_force,
                    'pre_load_displacement':d,
                    'energy_total': stored_energy,
                    'safety_factor': 1.5,
                    'E':E,
                    'A':A,
                    'I':I
                }

    return dict_flexure

if __name__ == "__main__":
    dict_flexure = parameters()
    print(dict_flexure)