#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Function to calculate the expectation value (mean) between two vectors
def expectation_value(vector1, vector2):
    dot_product = np.dot(vector1, vector2)
    expectation = dot_product / len(vector1)
    return expectation


#########################################################
# Generate Random Vectors (Remove for data)##############
# Create an array of vectors with the specified pattern
start = 5.0
end = 2.0
step_up = 0.0002
step_down = -0.0003

num_up_steps = int((7.0 - start) / step_up) + 1
num_down_steps = int((end - 7.0) / step_down) + 1

vector_length = 3
num_vectors = num_up_steps + num_down_steps

array_of_vectors = np.empty((num_vectors, vector_length))

# Fill the array with the desired pattern
current_value = start
for i in range(num_up_steps):
    array_of_vectors[i] = [current_value, current_value, current_value]
    current_value += step_up
    
for i in range(num_up_steps, num_vectors):
    array_of_vectors[i] = [current_value, current_value, current_value]
    current_value += step_down
        
vectors = array_of_vectors
############################################################

vectors = np.array(vectors)

# Vector 1 from the original data (Relaxed State)
vector1 = np.array([1, 1, 1])

# Calculate the expectation value between vector1 and each random vector
expectation_values = [expectation_value(vector1, vector) for vector in vectors]

# This will be our times
x_values = np.array(range(1, num_vectors + 1))

expectation_values = np.array(expectation_values)


###################################################
# Numerical Derivative

# Calculate the step size
step_size = x_values[1] - x_values[0]

# Calculate the numerical derivative of the expectation values
Phi = np.gradient(expectation_values, step_size)
####################################################


# Negative of the derivative of the fitted polynomial function
negative_Phi_derivative = -Phi

# Perform Fourier Transform on the negative derivative
# sampling rate indicates how frequently the signal is sampled or measured. The sampling rate 
#     is typically expressed in Hertz (Hz), which represents the number of samples per second.
sampling_rate = 2e+12  # calculated for every 0.5 ps
freqs = np.fft.fftfreq(len(negative_Phi_derivative), d=1/sampling_rate)
fft_values = np.fft.fft(negative_Phi_derivative)

# Get the real and imaginary parts separate
real_part = np.real(fft_values)
imaginary_part = np.imag(fft_values)


# Visualize this all


# Plot the original signal, its derivative, Fourier Transform, and the negative Phi derivative
plt.subplot(4, 1, 1)
plt.plot(x_values, expectation_values)
plt.title("Time Correlation Expectation")

plt.subplot(4, 1, 2)
plt.plot(x_values, negative_Phi_derivative)
plt.title("Negative Phi Derivative")

plt.subplot(4, 1, 3)
plt.plot(freqs, np.abs(fft_values))
plt.title("Fourier Transform of Negative Phi Derivative")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.xlim(-1e10, 1e10)

# Plot the real and imaginary parts of the Fourier transform
plt.figure()
plt.plot(freqs, real_part, label='Real Part')
plt.plot(freqs, imaginary_part, label='Imaginary Part')
plt.legend()
plt.title('Real and Imaginary Parts of Fourier Transform')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
# Set x-axis limit to start at 0
#plt.xlim(left=0, right=0.1)
plt.xlim(0, 0.1e11)

plt.tight_layout()
plt.show()

