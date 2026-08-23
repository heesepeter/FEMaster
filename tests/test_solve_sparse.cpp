#include "../src/solve/sparse/solve_sparse_direct.h"

#include <gtest/gtest.h>

using namespace fem;

TEST(SparseDirectSolver, SolvesSymmetricPositiveDefiniteMultipleRhs) {
    TripletList entries{
        {0, 0, 6.0}, {0, 1, 2.0}, {0, 2, 1.0},
        {1, 0, 2.0}, {1, 1, 5.0}, {1, 2, 2.0},
        {2, 0, 1.0}, {2, 1, 2.0}, {2, 2, 4.0}
    };
    SparseMatrix matrix(3, 3);
    matrix.setFromTriplets(entries.begin(), entries.end());

    DynamicMatrix expected(3, 2);
    expected << 1.0, -2.0,
                2.0,  0.5,
               -1.0,  3.0;
    const DynamicMatrix rhs = matrix * expected;

    const DynamicMatrix actual = solver::solve_direct(
        solver::CPU,
        matrix,
        rhs,
        solver::DirectSolverMatrixType::SPD);

    EXPECT_TRUE(actual.allFinite());
    EXPECT_TRUE(actual.isApprox(expected, 1e-12));
    EXPECT_LT((matrix * actual - rhs).norm(), 1e-11);
}
