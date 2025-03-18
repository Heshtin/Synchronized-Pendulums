
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>

// Function to compute the derivatives. Returns an array containing:
// [w1, w2, v, w1_2, w2_2, v_2]
std::array<double, 6> derv(double theta1, double theta2, double x,
                           double w1, double w2, double v,
                           double l1, double l2, double R,
                           double k, double c_x)
{
    const double g = 1.0; // Gravity constant

    double sin_theta1 = std::sin(theta1);
    double cos_theta1 = std::cos(theta1);
    double sin_theta2 = std::sin(theta2);
    double cos_theta2 = std::cos(theta2);

    double a = g * (sin_theta1 * cos_theta1 + sin_theta2 * cos_theta2);
    double b = l2 * (w2 * w2) * sin_theta2 + l1 * (w1 * w1) * sin_theta1;
    double c = R + 2 - (cos_theta1 * cos_theta1) - (cos_theta2 * cos_theta2);
    double v_2 = (a + b - c_x * v - k * x) / c; // dv/dt
    double w1_2 = (g * sin_theta1 + v_2 * cos_theta1) * (-1.0 / l1); // dw1/dt
    double w2_2 = (g * sin_theta2 + v_2 * cos_theta2) * (-1.0 / l2); // dw2/dt

    return {w1, w2, v, w1_2, w2_2, v_2};
}

// Runge-Kutta 4 integrator.
// Returns a vector of states. Each state is an array containing:
// [theta1, theta2, w1, w2, x, v, t]
std::vector<std::array<double, 7>> rk4(double theta1, double theta2, double w1,
                                       double w2, double x, double v,
                                       double l1, double l2, double R,
                                       double k, double c_x,
                                       double dt, int N)
{
    std::vector<std::array<double, 7>> result;
    result.reserve(N);

    for (int i = 0; i < N; i++) {
        // Save current state with time t = dt * i.
        result.push_back({theta1, theta2, w1, w2, x, v, dt * i});

        // Compute the RK4 slopes:
        auto f0 = derv(theta1, theta2, x, w1, w2, v, l1, l2, R, k, c_x);
        auto f1 = derv(theta1 + f0[0] * (dt / 2),
                       theta2 + f0[1] * (dt / 2),
                       x      + f0[2] * (dt / 2),
                       w1     + f0[3] * (dt / 2),
                       w2     + f0[4] * (dt / 2),
                       v      + f0[5] * (dt / 2),
                       l1, l2, R, k, c_x);
        auto f2 = derv(theta1 + f1[0] * (dt / 2),
                       theta2 + f1[1] * (dt / 2),
                       x      + f1[2] * (dt / 2),
                       w1     + f1[3] * (dt / 2),
                       w2     + f1[4] * (dt / 2),
                       v      + f1[5] * (dt / 2),
                       l1, l2, R, k, c_x);
        auto f3 = derv(theta1 + f2[0] * dt,
                       theta2 + f2[1] * dt,
                       x      + f2[2] * dt,
                       w1     + f2[3] * dt,
                       w2     + f2[4] * dt,
                       v      + f2[5] * dt,
                       l1, l2, R, k, c_x);

        // Update state using weighted average of slopes:
        theta1 += (f0[0] + 2 * f1[0] + 2 * f2[0] + f3[0]) * dt / 6;
        theta2 += (f0[1] + 2 * f1[1] + 2 * f2[1] + f3[1]) * dt / 6;
        x      += (f0[2] + 2 * f1[2] + 2 * f2[2] + f3[2]) * dt / 6;
        w1     += (f0[3] + 2 * f1[3] + 2 * f2[3] + f3[3]) * dt / 6;
        w2     += (f0[4] + 2 * f1[4] + 2 * f2[4] + f3[4]) * dt / 6;
        v      += (f0[5] + 2 * f1[5] + 2 * f2[5] + f3[5]) * dt / 6;
    }

    return result;
}

int main()
{
    // Example simulation using RK4

    // Initial conditions:
    double theta1 = 0.1;
    double theta2 = 0.7;
    double w1 = 0.0;  // Angular velocity of ball 1
    double w2 = 0.0;  // Angular velocity of ball 2
    double x = 0.0;   // Position of the center of the beam
    double v = 0.0;   // Velocity of the beam

    // System parameters:
    double l1 = 1.0;  // Length of pendulum 1
    double l2 = 1.0;  // Length of pendulum 2
    double R = 0.56;  // Big mass to ball mass ratio
    double k = 0.0;   // Spring constant
    double c_x = 0.0; // Damping constant

    double dt = 0.004;  // Time step
    int N = 400000;     // Number of steps

    // Run the integration:
    std::vector<std::array<double, 7>> simulationData = rk4(theta1, theta2, w1, w2, x, v,
                                                              l1, l2, R, k, c_x, dt, N);

    // (Optional) Print a few sample data points:
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "t, theta1, theta2, w1, w2, x, v\n";
    for (size_t i = 0; i < simulationData.size(); i += simulationData.size() / 10) {
        const auto& state = simulationData[i];
        std::cout << state[6] << ", " << state[0] << ", " << state[1] << ", "
                  << state[2] << ", " << state[3] << ", " << state[4] << ", "
                  << state[5] << "\n";
    }

    return 0;
}
