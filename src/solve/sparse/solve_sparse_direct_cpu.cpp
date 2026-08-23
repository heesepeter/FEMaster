#include "solve_sparse_direct.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../core/config.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

#ifdef USE_ACCELERATE_SPARSE
#include <Accelerate/Accelerate.h>
#include <vector>
#endif

#ifdef USE_MKL
#include <mkl.h>
#endif

namespace fem::solver::detail {

#ifdef USE_ACCELERATE_SPARSE
namespace {

bool solve_accelerate_spd(SparseMatrix& mat,
                          const DynamicMatrix& rhs,
                          DynamicMatrix& solution) {
    mat.makeCompressed();

    std::vector<long> column_starts(static_cast<std::size_t>(mat.cols()) + 1, 0);
    std::vector<int> row_indices;
    std::vector<double> values;
    row_indices.reserve(static_cast<std::size_t>(mat.nonZeros() / 2 + mat.rows()));
    values.reserve(row_indices.capacity());

    for (int column = 0; column < mat.outerSize(); ++column) {
        column_starts[static_cast<std::size_t>(column)] =
            static_cast<long>(row_indices.size());
        for (SparseMatrix::InnerIterator entry(mat, column); entry; ++entry) {
            if (entry.row() < column) continue;
            row_indices.push_back(entry.row());
            values.push_back(entry.value());
        }
    }
    column_starts.back() = static_cast<long>(row_indices.size());

    SparseAttributes_t attributes{};
    attributes.transpose = false;
    attributes.triangle = SparseLowerTriangle;
    attributes.kind = SparseSymmetric;

    SparseMatrixStructure structure{
        static_cast<int>(mat.rows()),
        static_cast<int>(mat.cols()),
        column_starts.data(),
        row_indices.data(),
        attributes,
        1
    };
    SparseMatrix_Double matrix{structure, values.data()};
    SparseOpaqueFactorization_Double factor =
        SparseFactor(SparseFactorizationCholesky, matrix);

    if (factor.status < SparseStatusOK) {
        SparseCleanup(factor);
        return false;
    }

    solution = rhs;
    DenseMatrix_Double dense{
        static_cast<int>(solution.rows()),
        static_cast<int>(solution.cols()),
        static_cast<int>(solution.outerStride()),
        SparseAttributes_t{},
        solution.data()
    };
    SparseSolve(factor, dense);
    SparseCleanup(factor);
    return solution.allFinite();
}

} // namespace
#endif

DynamicMatrix solve_direct_cpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type) {
    Timer t {};
    t.start();

    if (matrix_type == DirectSolverMatrixType::General) {
        DynamicMatrix sol;

#ifdef USE_MKL
        logging::info(true, "Using MKL PardisoLU solver");
        Eigen::PardisoLU<SparseMatrix> solver {};
        solver.compute(mat);
        if (solver.info() == Eigen::Success) {
            sol = solver.solve(rhs);
        }
#else
        logging::info(true, "Using Eigen SparseLU solver");
        Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver {};
        solver.compute(mat);
        if (solver.info() == Eigen::Success) {
            sol = solver.solve(rhs);
        }
#endif

        if (solver.info() != Eigen::Success) {
            logging::warning(false,
                "General sparse LU failed; falling back to SparseQR");
            Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
            qr.compute(mat);
            sol = qr.solve(rhs);
            logging::error(qr.info() == Eigen::Success,
                           "Solving general sparse system failed with SparseQR");
        }

        t.stop();
        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        logging::info(true, "residual    : ", (rhs - mat * sol).norm() / rhs.norm());
        return sol;
    }

#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
    int mkl_max_threads = mkl_get_max_threads();
    logging::info(true, "MKL max threads: ", mkl_max_threads);

    logging::info(true, "Using MKL PardisoLDLT solver");
    Eigen::PardisoLDLT<SparseMatrix> solver {};

    solver.compute(mat);
    logging::warning(solver.info() == Eigen::Success, "Decomposition failed with PardisoLDLT");
    DynamicMatrix sol = solver.solve(rhs);

    if (solver.info() != Eigen::Success) {
        logging::warning(true, "Solving failed with PardisoLDLT");
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
        qr.compute(mat);
        sol = qr.solve(rhs);
        logging::error(qr.info() == Eigen::Success, "Solving failed with SparseQR");
    }
#elif defined(USE_ACCELERATE_SPARSE)
    logging::info(true, "Using Apple Accelerate sparse Cholesky solver");
    DynamicMatrix sol;
    if (!solve_accelerate_spd(mat, rhs, sol)) {
        logging::warning(true,
            "Accelerate sparse Cholesky failed; falling back to Eigen SparseQR");
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
        qr.compute(mat);
        sol = qr.solve(rhs);
        logging::error(qr.info() == Eigen::Success,
                       "Solving failed with SparseQR");
    }
#else
    logging::info(true, "Using Eigen SimplicialLDLT solver");

    Eigen::SimplicialLDLT<SparseMatrix> solver {};
    solver.compute(mat);
    DynamicMatrix sol;
    if (solver.info() == Eigen::Success) {
        sol = solver.solve(rhs);
    }
    if (solver.info() != Eigen::Success) {
        logging::warning(true, "SimplicialLDLT failed; falling back to SparseQR");
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
        qr.compute(mat);
        sol = qr.solve(rhs);
        logging::error(qr.info() == Eigen::Success, "Solving failed with SparseQR");
    }
#endif

    t.stop();
    logging::info(true, "Solving finished");
    logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    logging::info(true, "residual    : ", (rhs - mat * sol).norm() / (rhs.norm()));

    return sol;
}

} // namespace fem::solver::detail
