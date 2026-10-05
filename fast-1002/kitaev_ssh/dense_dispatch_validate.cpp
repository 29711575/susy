#include "pure_projector_fast.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

double error(const MatType &left, const MatType &right) {
    return (left - right).cwiseAbs().maxCoeff();
}

} // namespace

int main() {
    try {
        const int dimension = 6;
        const int columns = 3;
        MatType dense(dimension, dimension);
        MatType q(dimension, columns);
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                dense(row, col) = DataType(
                    0.13 * double(1 + row + 2 * col),
                    -0.07 * double(1 + 2 * row + col));
            }
        }
        for (int col = 0; col < columns; ++col) {
            for (int row = 0; row < dimension; ++row) {
                q(row, col) = DataType(
                    0.11 * double(1 + row + col),
                    0.05 * double(1 + 2 * row - col));
            }
        }

        const PureFactorMatrix validatedDense =
            PureFactorMatrix::fromValidatedDense(dense);
        require(validatedDense.denseStorage() != nullptr,
                "validated dense factor did not retain dense storage");

        MatType actual;
        const MatType expected = dense * q;
        pureProjectorDensePropagationStats = {};
        pureProjectorPropagate(validatedDense, q, actual, false);
        require(error(actual, expected) == 0.0,
                "stable dense propagation changed matrix elements");
        pureProjectorPropagate(validatedDense, q, actual, true);
        require(error(actual, expected) == 0.0,
                "magnitude dense propagation changed matrix elements");
        require(pureProjectorDensePropagationStats.stable_direct_dense == 1,
                "stable dense dispatch count mismatch");
        require(pureProjectorDensePropagationStats.magnitude_direct_dense == 1,
                "magnitude dense dispatch count mismatch");
        require(pureProjectorDensePropagationStats.stable_expanded_factor == 0,
                "stable dense path expanded the factor");
        require(pureProjectorDensePropagationStats.magnitude_expanded_factor == 0,
                "magnitude dense path expanded the factor");

        MatType local = MatType::Identity(dimension, dimension);
        local(0, 1) = DataType(0.23, -0.11);
        local(1, 0) = DataType(-0.17, 0.09);
        const PureFactorMatrix generic(local);
        require(generic.denseStorage() == nullptr,
                "local generic factor unexpectedly retained dense storage");
        const MatType localExpected = local * q;
        pureProjectorDensePropagationStats = {};
        pureProjectorPropagate(generic, q, actual, false);
        require(error(actual, localExpected) == 0.0,
                "stable generic propagation changed matrix elements");
        pureProjectorPropagate(generic, q, actual, true);
        require(error(actual, localExpected) == 0.0,
                "magnitude generic propagation changed matrix elements");
        require(pureProjectorDensePropagationStats.stable_direct_dense == 0,
                "stable generic factor used dense dispatch");
        require(pureProjectorDensePropagationStats.magnitude_direct_dense == 0,
                "magnitude generic factor used dense dispatch");
        require(pureProjectorDensePropagationStats.stable_expanded_factor == 1,
                "stable generic dispatch count mismatch");
        require(pureProjectorDensePropagationStats.magnitude_expanded_factor == 1,
                "magnitude generic dispatch count mismatch");

        std::cout << "{\"stable_direct_dense\":"
                  << pureProjectorDensePropagationStats.stable_direct_dense
                  << ",\"stable_expanded_factor\":"
                  << pureProjectorDensePropagationStats.stable_expanded_factor
                  << ",\"magnitude_direct_dense\":"
                  << pureProjectorDensePropagationStats.magnitude_direct_dense
                  << ",\"magnitude_expanded_factor\":"
                  << pureProjectorDensePropagationStats.magnitude_expanded_factor
                  << ",\"status\":\"PASS\"}\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "dense_dispatch_validate FAIL: " << error.what() << '\n';
        return 1;
    }
}
