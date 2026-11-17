import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from parameters import parameters

def energy_at_deflection(theta, values):
    """Return axial, bending, and total strain energy at ``theta``."""
    flex_length_initial = values['flex_length_initial']
    gamma = values['gamma']
    E = values['E']
    A = values['A']
    I = values['I']
    K_Theta = values['K_Theta']

    L_prime_initial = flex_length_initial * gamma
    displacement_vertical=flex_length_initial * (1-np.cos(theta))
    flex_length_deformed = flex_length_initial -displacement_vertical
    L_diagonal=flex_length_deformed-(1-gamma)*flex_length_initial
    u=L_diagonal*np.sin(theta)
    L_prime_deformed = np.sqrt(
        u ** 2
        + (
            (gamma*2*flex_length_initial-displacement_vertical)
            / 2
        ) ** 2
    )
    axial_displacement = L_prime_deformed - L_prime_initial

    k_l = E * A / L_prime_initial
    axial_energy = 0.5 * k_l * axial_displacement ** 2
    Force_axial = (
        axial_energy / axial_displacement
        if axial_displacement != 0
        else 0.0
    )

    k_theta = 2 * gamma * K_Theta * E * I / (flex_length_initial)
    bending_energy = 2 * 0.5 * k_theta * theta ** 2
    Force_bending = bending_energy / theta if theta != 0 else 0.0
    Force_total =   Force_bending + values['preload_force']

    total_energy = bending_energy

    #print(Force_axial, Force_bending, Force_total, axial_energy,axial_displacement,L_prime_initial,k_l)

    return axial_energy, bending_energy, total_energy, Force_total


def find_deflection_angle(values, theta_limit=np.radians(50), samples=1000):
    """Solve total strain energy(theta) = required energy."""
    target_energy = values['energy_total']
    initial_angle = values['init_angle']
    theta_grid = np.linspace(0.0, theta_limit-initial_angle, samples + 1)
    energy_grid = np.array(
        [energy_at_deflection(theta, values)[2] for theta in theta_grid]
    )
    residual_grid = energy_grid - target_energy

    for index in range(samples):
        if residual_grid[index] == 0:
            theta = theta_grid[index]
            break
        if residual_grid[index] * residual_grid[index + 1] >0:
            continue

        lower = theta_grid[index]
        upper = theta_grid[index + 1]
        for _ in range(60):
            midpoint = (lower + upper) / 2
            midpoint_residual = (
                energy_at_deflection(midpoint, values)[2] - target_energy
            )
            if residual_grid[index] * midpoint_residual <= 0:
                upper = midpoint
            else:
                lower = midpoint
        theta = (lower + upper) / 2
        break
    else:
        return None

    axial_energy, bending_energy, total_energy, Force_total = energy_at_deflection(
        theta, values
    )
    return {
        'deflection_angle': theta,
        'axial_energy': axial_energy,
        'bending_energy': bending_energy,
        'total_energy': total_energy,
        'energy_residual': total_energy - target_energy,
        'Force_total': Force_total,
    }


def find_perfect_combination(parameters_dict):
    best_combination = None
    smallest_residual = float('inf')

    for combination, values in parameters_dict.items():
        result = find_deflection_angle(values)
        if result is None:
            continue
        residual = abs(result['energy_residual'])
        if residual < smallest_residual:
            smallest_residual = residual
            best_combination = {
                'combination': combination,
                'preload_force': values['preload_force'],
                **result,
            }
    return best_combination

def dtheta_dy(theta, vertical_displacement):
    """Return the local slope using the supplied angle and displacement units."""
    return np.gradient(theta, vertical_displacement)


def plot_deflection_angle_over_time(
    parameters_dict, duration=2 * np.pi, samples=100
):
    """Plot angle during linear compression of the best combination."""
    best_combination = find_perfect_combination(parameters_dict)
    if best_combination is None:
        raise ValueError(
            'No parameter combination can store the required energy '
            'within the deflection limit.'
        )

    values = parameters_dict[best_combination['combination']]
    flex_length = values['flex_length_initial']
    initial_angle = values['init_angle']
    preload_displacement = abs(values['pre_load_displacement'])
    maximum_angle = best_combination['deflection_angle']
    maximum_vertical_displacement = flex_length * (1 - np.cos(maximum_angle))

    time = np.linspace(0.0, duration, samples)
    compression_displacement = np.linspace(
        0.0, maximum_vertical_displacement, samples
    )
    compression_angle = np.arccos(
        1 - compression_displacement / flex_length
    )
    deflection_angle = initial_angle + compression_angle
    deflection_angle_degrees = np.degrees(deflection_angle)
    vertical_displacement = (
        preload_displacement + compression_displacement
    )
    vertical_displacement_mm = vertical_displacement * 1000
    flex_length_deformed = flex_length - compression_displacement
    diagonal_length = (
        flex_length_deformed
        - (1 - values['gamma']) * flex_length
    )
    horizontal_displacement = diagonal_length * np.sin(compression_angle)
    horizontal_displacement_mm = horizontal_displacement * 1000
    total_energy = np.array(
        [
            energy_at_deflection(angle, values)[2]
            for angle in compression_angle
        ]
    )
    Force_tot = np.array(
        [
            energy_at_deflection(angle, values)[3]
            for angle in compression_angle
        ]
    )
    dtheta_dy_over_time = np.array(
        dtheta_dy(deflection_angle_degrees, vertical_displacement_mm)
    )

    figure, axes = plt.subplots(2, 3, figsize=(16, 8))
    axes[0, 0].plot(time, deflection_angle_degrees, linewidth=2)
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Deflection angle (degrees)')
    axes[0, 0].set_title('Deflection angle')
    axes[0, 0].grid(True)

    axes[0, 1].plot(time, total_energy, linewidth=2)
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('Total energy (J)')
    axes[0, 1].set_title('Total stored energy')
    axes[0, 1].grid(True)

    axes[1, 0].plot(time, vertical_displacement_mm, linewidth=2)
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Vertical displacement (mm)')
    axes[1, 0].set_title('Vertical compression')
    axes[1, 0].grid(True)

    axes[1, 1].plot(time, horizontal_displacement_mm, linewidth=2)
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Horizontal displacement (mm)')
    axes[1, 1].set_title('Horizontal displacement, u')
    axes[1, 1].grid(True)

    axes[0, 2].plot(time, Force_tot, linewidth=2)
    axes[0, 2].set_xlabel('Time (s)')
    axes[0, 2].set_ylabel('Force (N)')
    axes[0, 2].set_title('Total force')
    axes[0, 2].grid(True)

    axes[1, 2].plot(
        vertical_displacement_mm, deflection_angle_degrees, linewidth=2
    )
    axes[1, 2].set_xlim(0.0, vertical_displacement_mm[-1])
    axes[1, 2].set_xlabel('Vertical displacement (mm)')
    axes[1, 2].set_ylabel('Deflection angle (degrees)')
    axes[1, 2].set_title('Deflection angle over vertical displacement')
    axes[1, 2].grid(True)

    figure.suptitle(f"Compression response: {best_combination['combination']}")
    figure.tight_layout()
    plt.show()

    best_combination['time'] = time
    best_combination['vertical_displacement'] = vertical_displacement
    best_combination['horizontal_displacement'] = horizontal_displacement
    best_combination['deflection_angle_over_time'] = deflection_angle
    best_combination['total_energy_over_time'] = total_energy
    best_combination['Force_tot'] = Force_tot
    best_combination['dtheta_dy_over_time'] = dtheta_dy_over_time
    return best_combination


if __name__ == '__main__':
    result = plot_deflection_angle_over_time(parameters())
    result_rows = [
        ['Flexure Width', f"{(result['combination'][0]*1000):.2f} mm"],
        ['Flexure Thickness', f"{(result['combination'][1]*1000):.2f} mm"],
        ['Starting Deflection angle', f"{result['combination'][2]:.2f} degrees"],
        ['Deflection angle', f"{np.degrees(result['deflection_angle']):.2f} degrees"],
        ['Pre load force', f"{result['preload_force']:.2f} N"],
        ['Total force', f"{result['Force_total']:.2f} N"],
        [
            'Applied force',
            f"{result['Force_total'] - result['preload_force']:.2f} N",
        ],
        [
                    'Total Applied force per Leg',
                    f"{(result['Force_total'] - result['preload_force'])*4:.2f} N",
                ],
        [
            'Total vertical deflection',
            f"{2 * np.max(result['vertical_displacement']) * 1000:.2f} mm",
        ],
        [
            'Average dtheta/dy',
            f"{np.mean(result['dtheta_dy_over_time']):.2f} deg/mm",
        ],
    ]
    red = '\033[1;91m'
    reset = '\033[0m'
    print(
        f"\n{red}BEST COMBINATION{reset}\n"
        + tabulate(
            result_rows,
            headers=['Quantity', 'Value'],
            tablefmt='fancy_grid',
        )
    )