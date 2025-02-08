
import numpy as np
from scipy.odr import ODR, Model, RealData


def load_data_with_commas(filename):
    # Read the file, replace commas with dots, and load as a NumPy array
    with open(filename, 'r') as file:
        data = file.read().replace(',', '.')
    return np.loadtxt(data.splitlines())

# Load the data
data = np.loadtxt("C:\Users\Ale\Desktop\f_B.txt")  # Ensure the file is formatted as: x, err(x), y, err(y)

# Extract columns
x = data[:, 0]         # x
err_x = data[:, 1]     # Error in x
y = data[:, 2]         # y
err_y = data[:, 3]     # Error in y

# Define the linear model
def linear_model(params, x):
    m, c = params  # m: slope, c: intercept
    return m * x + c

# Prepare the data for ODR
model = Model(linear_model)
data = RealData(x, y, sx=err_x, sy=err_y)
odr = ODR(data, model, beta0=[1, 0])  # Initial guess for m and c

# Run the fit
output = odr.run()

# Extract fit parameters
m, c = output.beta  # Slope and intercept
print(f"Slope (m): {m}")
print(f"Intercept (c): {c}")

# Extract standard errors
m_err, c_err = output.sd_beta
print(f"Error in Slope (m): {m_err}")
print(f"Error in Intercept (c): {c_err}")
