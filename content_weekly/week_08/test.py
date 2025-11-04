# Import required libraries
# NOTE: Nothing for you to do here
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import warnings
warnings.filterwarnings('ignore')

pi = 3.141592653589793

def velocity_source_sink(Q, x0, y0, x, y):
    """
    Calculate velocity components for a source/sink
    
    Parameters:
    Q: Source/sink strength (positive for source, negative for sink)
    x0, y0: Source/sink location
    X, Y: Coordinate grids
    
    Returns:
    u, v: Velocity components
    """
    r_squared = (x - x0)**2 + (y - y0)**2
    # Avoid division by zero, don't change this line
    r_squared = np.maximum(r_squared, 1e-10)
    
    u = (Q * (x - x0)) / (2 * pi * r_squared)
    v = (Q * (y - y0)) / (2 * pi * r_squared)
    
    return u, v

def stream_function_source_sink(Q, x0, y0, x, y):
    """
    Calculate stream function for a source/sink

    Parameters:
    Q: Source/sink strength (positive for source, negative for sink)
    x0, y0: Source/sink location
    X, Y: Coordinate grids

    Returns:
    psi: Stream function
    """
    theta = np.arctan2(y - y0, x - x0) # Note, use np.arctan2 for correct angle calculation
    psi = (Q / (2 * pi)) * theta
    return psi

def velocity_potential_source_sink(Q, x0, y0, x, y):
    """
    Calculate velocity potential for a source/sink

    Parameters:
    Q: Source/sink strength (positive for source, negative for sink)
    x0, y0: Source/sink location
    X, Y: Coordinate grids
    
    Returns:
    phi: Velocity potential
    """
    r = ((x - x0)**2 + (y - y0)**2)** 0.5
    # Avoid log(0), don't change this line
    r = np.maximum(r, 1e-10)
    phi = Q / (2 * pi) * np.log(r)
    return phi

# Define the computational domain
x_min, x_max = -3.0, 3.0
y_min, y_max = -2.0, 2.0
nx, ny = 200, 150

# Define source and sink parameters
Q_source = 5
Q_sink = -5
x_source, y_source = (-1, 0)
x_sink, y_sink = (1, 0)

# Create coordinate grids
x = np.linspace(x_min, x_max, nx)
y = np.linspace(y_min, y_max, ny)
x, y = np.meshgrid(x, y)

# Calculate velocity components
u_source, v_source = velocity_source_sink(Q_source, x_source, y_source, x, y)
u_sink, v_sink = velocity_source_sink(Q_sink, x_sink, y_sink, x, y)

# Superimpose the velocity fields
u_total = u_source + u_sink
v_total = v_source + v_sink
velocity_magnitude = (u_total**2 + v_total**2)**0.5

# Create a mask to avoid singularities near source/sink
# NOTE: Nothing to change here, this will make the visualization cleaner
mask = ((x - x_source)**2 + (y - y_source)**2 < 0.01) | ((x - x_sink)**2 + (y - y_sink)**2 < 0.01)
velocity_magnitude = np.ma.masked_where(mask, velocity_magnitude)

# Calculate stream function and velocity potential
psi_source = stream_function_source_sink(Q_source, x_source, y_source, x, y)
psi_sink = stream_function_source_sink(Q_sink, x_sink, y_sink, x, y)
psi_total = psi_source + psi_sink

phi_source = velocity_potential_source_sink(Q_source, x_source, y_source, x, y)
phi_sink = velocity_potential_source_sink(Q_sink, x_sink, y_sink, x, y)
phi_total = phi_source + phi_sink
# Create comprehensive flow net visualization
# Use the plt.subplots method to create figure and 2 axes
# TODO: Your code here for creating figure and axes

fig, axes = plt.subplots(2, 1, figsize=(15, 12))
fig.suptitle('Groundwater Flow Around a Source and Sink', fontsize=16, fontweight='bold')

# Stream function contours/flow net
# Use the ax.streamplot, ax.contour, and ax.scatter methods to plot
# the streamlines, equipotential lines, and mark the source/sink locations respectively
# TODO: Your code here for plotting

ax1 = axes[0]

ax1.scatter(x_source, y_source, s=100, c='black', marker='o', label='Source', zorder = 3)
ax1.scatter(x_sink, y_sink, s=100, c='black', marker='s', label='Sink', zorder = 3)

ax1.streamplot(x, y, u_total, v_total, linewidth = 1, color = 'dodgerblue', zorder = 1)

cs = ax1.contour(x, y, phi_total, colors='red', levels=20, zorder = 2)
ax1.clabel(cs, fontsize=10, fmt=lambda x: f"{x:.2f} m", zorder = 2)

ax1.legend()
ax1.set_title("Groundwater Flow Field")
ax1.set_xlabel("X Coordinate (m)")
ax1.set_ylabel("Y Coordinate (m)")

# Velocity magnitude
# Use the ax.contourf and ax.scatter methods to plot
# the velocity magnitude and mark the source/sink locations respectively
# TODO: Your code here for plotting

ax2 = axes[1]

ax2.scatter(x_source, y_source, s=100, c='black', marker='o', label='Source', zorder = 3)
ax2.scatter(x_sink, y_sink, s=100, c='black', marker='s', label='Sink', zorder = 3)

cs = ax2.contourf(x, y, velocity_magnitude, cmap='GnBu', levels=20, zorder = 2)

ax2.legend()
cbar = plt.colorbar(cs, ax=ax2)
cbar.set_label('Groundwater Velocity (m/s)')
ax2.set_title("Groundwater Velocity Distribution")
ax2.set_xlabel("X Coordinate (m)")
ax2.set_ylabel("Y Coordinate (m)")

# Finalize and show plots, tight layout fits everything nicely
plt.tight_layout()

def flow_around_cylinder(U, R, r, theta):
    """
    Calculate velocity components for flow around a cylinder
    
    Parameters:
    U: Free stream velocity
    R: Cylinder radius
    r: Radial distance from cylinder center
    theta: Angular coordinate
    
    Returns:
    ur, utheta: Radial and tangential velocity components
    """
    # Avoid division by zero and points inside cylinder
    # NOTE: Don't change this line
    r_safe = np.maximum(r, R + 1e-6)
    
    # Radial velocity component
    ur = U * np.cos(theta) * (1 - (R**2)/(r**2))
    # Tangential velocity component
    utheta = -U * np.sin(theta) * (1 - (R**2)/(r**2))
    return ur, utheta

def stream_function_cylinder(U, R, r, theta):
    """
    Calculate stream function for flow around a cylinder
    """
    # NOTE: Don't change this line
    r_safe = np.maximum(r, R + 1e-6)
    psi = U * r * np.sin(theta) - (U * R**2 * np.sin(theta) / r)
    return psi

def velocity_potential_cylinder(U, R, r, theta):
    """
    Calculate velocity potential for flow around a cylinder
    """
    # NOTE: Don't change this line
    r_safe = np.maximum(r, R + 1e-6)
    phi = U * r * np.cos(theta) - (U * R**2 * np.cos(theta) / r)
    return phi


# Flow around a cylinder
# Define parameters
U = 2.0  # Free stream velocity
R = 1.0      # Cylinder radius
x_center, y_center = 0.0, 0.0  # Cylinder center

# Create coordinate grids
x_cyl = np.linspace(-R-0.5, R+0.5, 100)
y_cyl = np.linspace(-R-0.5, R+0.5, 100)
X_cyl, Y_cyl = np.meshgrid(x_cyl, y_cyl)

# Convert to polar coordinates relative to cylinder center
# NOTE: I'm giving this transoformation to you, don't change the next two lines
r = np.sqrt((X_cyl - x_center)**2 + (Y_cyl - y_center)**2)
theta = np.arctan2(Y_cyl - y_center, X_cyl - x_center)

# Calculate velocity components
ur, utheta = flow_around_cylinder(U, R, r, theta)

# Convert to Cartesian coordinates
# NOTE: I'm giving this transoformation to you, don't change the next two lines
u_cyl = ur * np.cos(theta) - utheta * np.sin(theta)
v_cyl = ur * np.sin(theta) + utheta * np.cos(theta)

# Calculate stream function and velocity potential
psi_cyl = stream_function_cylinder(U, R, r, theta)
phi_cyl = velocity_potential_cylinder(U, R, r, theta)

# Mask points inside the cylinder
# NOTE: Nothing to change here, this will make the visualization cleaner
mask_inside = r < R
u_cyl = np.ma.masked_where(mask_inside, u_cyl)
v_cyl = np.ma.masked_where(mask_inside, v_cyl)
psi_cyl = np.ma.masked_where(mask_inside, psi_cyl)
phi_cyl = np.ma.masked_where(mask_inside, phi_cyl)

# Pressure distribution (using Bernoulli's equation)
# TODO: Your code here for calculating
velocity_magnitude_cyl = (u_cyl**2 + v_cyl**2)**0.5
rho = 1
p_inf = 0
p = p_inf + 0.5 * rho * (U**2 - velocity_magnitude_cyl**2)


# Pressure coefficient Cp = 1 - (V/U)²
# TODO: Your code here for calculating
Cp = 1 - (velocity_magnitude_cyl / U)**2
# NOTE: Nothing to change here, this will make the visualization cleaner
Cp = np.ma.masked_where(mask_inside, Cp)
# Visualize flow around cylinder
# Create figure and 2 axes using plt.subplots
# TODO: Your code here for creating figure and axes
fig, axes = plt.subplots(2, 1, figsize=(15, 12))
fig.suptitle('Groundwater Flow Around a Cylinder', fontsize=16, fontweight='bold')

# Combined flow net
# Use the ax.streamplot, ax.contour to plot
# TODO: Your code here for plotting

ax1 = axes[0]

ax1.streamplot(x_cyl, y_cyl, u_cyl, v_cyl, linewidth = 1, color = 'dodgerblue', zorder = 1)

cs = ax1.contour(x_cyl, y_cyl, phi_cyl, colors='red', levels=10, zorder = 2)
ax1.clabel(cs, fontsize=10, fmt=lambda x: f"{x:.2f} m", zorder = 2)

ax1.legend()
ax1.set_title("Groundwater Flow Field")
ax1.set_xlabel("X Coordinate (m)")
ax1.set_ylabel("Y Coordinate (m)")

theta_circle = np.linspace(0, 2*np.pi, 200)
x_circle = x_center + R * np.cos(theta_circle)
y_circle = y_center + R * np.sin(theta_circle)
ax1.plot(x_circle, y_circle, color='gray', linewidth=2, zorder=3)

# Pressure coefficient distribution
# Use the ax.contourf to plot, and add colorbar
# TODO: Your code here for plotting

ax2 = axes[1]

cs = ax2.contourf(x_cyl, y_cyl, Cp, cmap='GnBu', levels=20, zorder = 2)

theta_circle = np.linspace(0, 2*np.pi, 200)
x_circle = x_center + R * np.cos(theta_circle)
y_circle = y_center + R * np.sin(theta_circle)
ax2.plot(x_circle, y_circle, color='gray', linewidth=2, zorder=3)

ax2.legend()
cbar = plt.colorbar(cs, ax=ax2)
cbar.set_label('Presure Coefficient')
ax2.set_title("Pressure Coefficient Distribution")
ax2.set_xlabel("X Coordinate (m)")
ax2.set_ylabel("Y Coordinate (m)")

def kirkham_head(D=1.0, B=1.0, D_c=0.2, Nxz=101, Nn=201):
    """
    Compute the steady-state hydraulic head H(x,z) in a rectangular aquifer
    using the analytical solution of Kirkham 

    Parameters
    ----------
    D : float
        Depth of the aquifer.
    B : float
        Breadth (horizontal extent) of the aquifer.
    D_c : float
        Water level at the open ditch.
    Nxz : int
        Number of discretization steps for x and z.
    Nn : int
        Number of terms to include in the series expansion (odd terms only).

    Returns
    -------
    H : ndarray of shape (Nxz, Nxz)
        Steady-state head field, with H[i, j] corresponding to (x[i], z[j]).
    x : ndarray
        1D array of x-coordinates (horizontal direction).
    z : ndarray
        1D array of z-coordinates (vertical direction).
    """
    # Discretization grids
    x = np.linspace(0,B,Nxz)
    z = np.linspace(0,D,Nxz)

    # Set up meshgrid for calculations
    # NOTE: Don't change this line
    x, z = np.meshgrid(x, z, indexing="ij")

    # Odd series terms
    # Use np.arange, and don't forget about the odd terms only!
    n_values = np.arange(1, Nn, 2)
    C = (8 * D) / math.pi**2

    # Initialize head field with zeros, of size (Nxz, Nxz)
    H = np.zeros((Nxz, Nxz))

    # Compute series sum term by term
    # TODO Your code here
    for n in n_values:
        lambda_n = n * math.pi / (2 * B)
        term = (1 / n ** 2) * np.sin(lambda_n * x) * (np.cosh(lambda_n * z) / np.cosh(lambda_n * D))
        H += term

    # NOTE: Nothing to do after this line
    H = D - H
    return H, x, z
# Kirkham solution - Define parameters
D = 10.0      # Aquifer depth (m)
B = 50.0      # Aquifer breadth (m)
Nxz = 151     # Grid resolution
Nn = 50

# Calculate head distributions for two different D_c values
H1, x_kirk, z_kirk = kirkham_head(D, B, 2.0, Nxz, Nn)
H2, _, _ = kirkham_head(D, B, 6.0, Nxz, Nn)

# Calculate the difference between the two scenarios
H_diff = H2 - H1

print("we're here")

print(x_kirk.flatten())

# Create meshgrids for plotting
# NOTE: Don't change these lines
X_kirk, Z_kirk = np.meshgrid(x_kirk.flatten(), z_kirk.flatten(), indexing='ij')

print(X_kirk[:10,:10])