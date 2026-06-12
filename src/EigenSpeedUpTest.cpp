// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppEigen.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <string>
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace {
  using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
  using Vec   = Eigen::VectorXd;

  inline void write_target_result(
      double* out_ptr,
      int out_nrow,
      int n,
      int m_targets,
      int jj,          // target position in `targets`
      int t,           // 0-based state index
      const Vec& h,
      const Vec& z,
      bool targets_in_rows,
      double hit_tol
  ) {
    for (int i = 0; i < n; ++i) {
      double val;
      if (i == t) {
        val = 0.0;
      } else {
        const double hi = h[i];
        const double zi = z[i];
        val = (!std::isfinite(hi) || !std::isfinite(zi) || std::abs(hi) <= hit_tol)
          ? NA_REAL
          : (hi * zi);   // exact legacy column-sum quantity
      }

      if (targets_in_rows) {
        // out is (m_targets x n): row = target, col = state
        out_ptr[jj + out_nrow * i] = val;
      } else {
        // out is (n x m_targets): row = state, col = target
        out_ptr[i + out_nrow * jj] = val;
      }
    }
  }

  inline std::string compute_one_target_flighttime(
      int t,
      int jj,
      SpMat& A,
      const std::vector<std::vector<int>>& row_slots,
      const std::vector<int>& diag_pos,
      const std::vector<double>& base_values,
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>>& solver,
      double* out_ptr,
      int out_nrow,
      int n,
      int m_targets,
      bool targets_in_rows,
      double hit_tol
  ) {
    double* Ax = A.valuePtr();

    // Impose row t = e_t^T without changing sparsity pattern.
    for (int pos : row_slots[t]) Ax[pos] = 0.0;
    Ax[diag_pos[t]] = 1.0;

    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
      for (int pos : row_slots[t]) Ax[pos] = base_values[pos];
      return std::string("SparseLU factorize failed at target ") +
        std::to_string(t + 1) + ": " + solver.lastErrorMessage();
    }

    // h = M_t^{-1} e_t
    Vec rhs1 = Vec::Zero(n);
    rhs1[t] = 1.0;
    Vec h = solver.solve(rhs1);

    // z = M_t^{-T} d, where d_i = 1/h_i for i != t, d_t = 0
    Vec rhs2 = Vec::Zero(n);
    for (int i = 0; i < n; ++i) {
      if (i != t && std::isfinite(h[i]) && std::abs(h[i]) > hit_tol) {
        rhs2[i] = 1.0 / h[i];
      }
    }
    Vec z = solver.transpose().solve(rhs2);

    write_target_result(
      out_ptr, out_nrow, n, m_targets, jj, t, h, z,
      targets_in_rows, hit_tol
    );

    // Restore original numeric values for next target.
    for (int pos : row_slots[t]) Ax[pos] = base_values[pos];

    return "";
  }
}

//' Eigen Matrix Multiplication
//'
//' An alternative to \code{format.POSIXct} based on the CCTZ library. The
//' \code{formatDouble} variant uses two vectors for seconds since the epoch
//' and fractional nanoseconds, respectively, to provide fuller resolution.
//'
//' @title Format a Datetime vector as a string vector
//' @param dtv A Datetime vector object to be formatted
//' @param fmt A string with the format, which is based on \code{strftime} with some
//'   extensions; see the CCTZ documentation for details.
//' @param lcltzstr The local timezone object for creation the CCTZ timepoint
//' @param tgttzstr The target timezone for the desired format
//' @return A string vector with the requested format of the datetime objects
//' @section Note:
//' Windows is now supported via the \code{g++-4.9} compiler, but note
//' that it provides an \emph{incomplete} C++11 library. This means we had
//' to port a time parsing routine, and that string formatting is more
//' limited. As one example, CCTZ frequently uses \code{"\%F \%T"} which do
//' not work on Windows; one has to use \code{"\%Y-\%m-\%d \%H:\%M:\%S"}.
//' @author Dirk Eddelbuettel
//' @examples
//' \dontrun{
//' now <- Sys.time()
//' formatDatetime(now)            # current (UTC) time, in full precision RFC3339
//' formatDatetime(now, tgttzstr="America/New_York")  # same but in NY
//' formatDatetime(now + 0:4)    # vectorised
//' }
// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}
//' Eigen Matrix Multiplication Faster
//'
//' An alternative to \code{format.POSIXct} based on the CCTZ library. The
//' \code{formatDouble} variant uses two vectors for seconds since the epoch
//' and fractional nanoseconds, respectively, to provide fuller resolution.
//'
//' @title Format a Datetime vector as a string vector
//' @param dtv A Datetime vector object to be formatted
//' @param fmt A string with the format, which is based on \code{strftime} with some
//'   extensions; see the CCTZ documentation for details.
//' @param lcltzstr The local timezone object for creation the CCTZ timepoint
//' @param tgttzstr The target timezone for the desired format
//' @return A string vector with the requested format of the datetime objects
//' @section Note:
//' Windows is now supported via the \code{g++-4.9} compiler, but note
//' that it provides an \emph{incomplete} C++11 library. This means we had
//' to port a time parsing routine, and that string formatting is more
//' limited. As one example, CCTZ frequently uses \code{"\%F \%T"} which do
//' not work on Windows; one has to use \code{"\%Y-\%m-\%d \%H:\%M:\%S"}.
//' @author Michael Bowman
//' @examples
//' \dontrun{
//' now <- Sys.time()
//' formatDatetime(now)            # current (UTC) time, in full precision RFC3339
//' formatDatetime(now, tgttzstr="America/New_York")  # same but in NY
//' formatDatetime(now + 0:4)    # vectorised
//' }
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}
//' Eigen Matrix Inverse
//'
//' An alternative to \code{format.POSIXct} based on the CCTZ library. The
//' \code{formatDouble} variant uses two vectors for seconds since the epoch
//' and fractional nanoseconds, respectively, to provide fuller resolution.
//'
//' @title Format a Datetime vector as a string vector
//' @param dtv A Datetime vector object to be formatted
//' @param fmt A string with the format, which is based on \code{strftime} with some
//'   extensions; see the CCTZ documentation for details.
//' @param lcltzstr The local timezone object for creation the CCTZ timepoint
//' @param tgttzstr The target timezone for the desired format
//' @return A string vector with the requested format of the datetime objects
//' @section Note:
//' Windows is now supported via the \code{g++-4.9} compiler, but note
//' that it provides an \emph{incomplete} C++11 library. This means we had
//' to port a time parsing routine, and that string formatting is more
//' limited. As one example, CCTZ frequently uses \code{"\%F \%T"} which do
//' not work on Windows; one has to use \code{"\%Y-\%m-\%d \%H:\%M:\%S"}.
//' @author Michael Bowman
//' @examples
//' \dontrun{
//' now <- Sys.time()
//' formatDatetime(now)            # current (UTC) time, in full precision RFC3339
//' formatDatetime(now, tgttzstr="America/New_York")  # same but in NY
//' formatDatetime(now + 0:4)    # vectorised
//' }
// [[Rcpp::export]]
SEXP eigenMatInverse(const Eigen::Map<Eigen::MatrixXd> A){
//Eigen::setNbThreads(8);
 Eigen::MatrixXd C = A.inverse();
 return Rcpp::wrap(C);
}
//' Eigen Matrix Inverse and Double Multiply
//'
//' An alternative to \code{format.POSIXct} based on the CCTZ library. The
//' \code{formatDouble} variant uses two vectors for seconds since the epoch
//' and fractional nanoseconds, respectively, to provide fuller resolution.
//'
//' @title Format a Datetime vector as a string vector
//' @param dtv A Datetime vector object to be formatted
//' @param fmt A string with the format, which is based on \code{strftime} with some
//'   extensions; see the CCTZ documentation for details.
//' @param lcltzstr The local timezone object for creation the CCTZ timepoint
//' @param tgttzstr The target timezone for the desired format
//' @return A string vector with the requested format of the datetime objects
//' @section Note:
//' Windows is now supported via the \code{g++-4.9} compiler, but note
//' that it provides an \emph{incomplete} C++11 library. This means we had
//' to port a time parsing routine, and that string formatting is more
//' limited. As one example, CCTZ frequently uses \code{"\%F \%T"} which do
//' not work on Windows; one has to use \code{"\%Y-\%m-\%d \%H:\%M:\%S"}.
//' @author Michael Bowman
//' @examples
//' \dontrun{
//' now <- Sys.time()
//' formatDatetime(now)            # current (UTC) time, in full precision RFC3339
//' formatDatetime(now, tgttzstr="America/New_York")  # same but in NY
//' formatDatetime(now + 0:4)    # vectorised
//' }
// [[Rcpp::export]]
SEXP eigenMatXMatXMatINV(const Eigen::Map<Eigen::MatrixXd> A,const Eigen::Map<Eigen::MatrixXd> B,const Eigen::Map<Eigen::MatrixXd> C){
  Eigen::MatrixXd Output = (A.inverse())*B*C;
  return Rcpp::wrap(Output);

}
//' Eigen for OutputFlightTime
//'
//' An alternative to \code{format.POSIXct} based on the CCTZ library. The
//' \code{formatDouble} variant uses two vectors for seconds since the epoch
//' and fractional nanoseconds, respectively, to provide fuller resolution.
//'
//' @title Format a Datetime vector as a string vector
//' @param dtv A Datetime vector object to be formatted
//' @param fmt A string with the format, which is based on \code{strftime} with some
//'   extensions; see the CCTZ documentation for details.
//' @param lcltzstr The local timezone object for creation the CCTZ timepoint
//' @param tgttzstr The target timezone for the desired format
//' @return A string vector with the requested format of the datetime objects
//' @section Note:
//' Windows is now supported via the \code{g++-4.9} compiler, but note
//' that it provides an \emph{incomplete} C++11 library. This means we had
//' to port a time parsing routine, and that string formatting is more
//' limited. As one example, CCTZ frequently uses \code{"\%F \%T"} which do
//' not work on Windows; one has to use \code{"\%Y-\%m-\%d \%H:\%M:\%S"}.
//' @author Michael Bowman
//' @examples
//' \dontrun{
//' now <- Sys.time()
//' formatDatetime(now)            # current (UTC) time, in full precision RFC3339
//' formatDatetime(now, tgttzstr="America/New_York")  # same but in NY
//' formatDatetime(now + 0:4)    # vectorised
//' }
// [[Rcpp::export]]
SEXP eigenOutputFlightTime(const Eigen::Map<Eigen::MatrixXd> Qmat,const Eigen::Map<Eigen::MatrixXd> Rmat){
  Eigen::MatrixXd N = (Eigen::MatrixXd::Identity(Qmat.rows(),Qmat.cols())-Qmat).inverse();
  Eigen::VectorXd B = N*Rmat;
  //Rcpp::Rcout<<(B.asDiagonal().toDenseMatrix()).inverse()<<std::endl;
  Eigen::MatrixXd temp_output = ((B.asDiagonal().toDenseMatrix()).inverse())*N*(B.asDiagonal().toDenseMatrix());
  // People and literature is incorrect. They are doing outgoing, NOT incoming.
  // We need to sum on columns to preserve incoming time. 
  // This means CMFPT going from s_t into s_t+1 (columns).
  return Rcpp::wrap(temp_output.colwise().sum());
 
  // NOT THIS becuase this is just saying how much would a state have outgoing
  //return Rcpp::wrap(temp_output.rowwise().sum());
}
/* this is column norm for fasr
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const int n = Q.rows();
  MatrixXd A = MatrixXd::Identity(n, n);
  A.noalias() -= Q;

  Eigen::PartialPivLU<MatrixXd> luA(A);
  Eigen::PartialPivLU<MatrixXd> luAT(A.transpose());

  VectorXd h = luA.solve(r);
  VectorXd y = luAT.solve(h.cwiseInverse());

  VectorXd out = (h.array() * y.array()).matrix();

  const double eps = 1e-15;
  for (int i = 0; i < n; ++i) {
    if (std::abs(h[i]) < eps) out[i] = NA_REAL;
  }

  return out;
*/

//Rcpp::NumericMatrix cmfpt_to

// [[Rcpp::export]]
Rcpp::NumericMatrix flighttime_sparse_targets_cpp(
    const Eigen::MappedSparseMatrix<double>& P,
    Rcpp::IntegerVector targets,
    bool targets_in_rows = true,
    double hit_tol = 1e-14,
    int n_threads = 1
) {
  const int n = P.rows();
  if (P.cols() != n) {
    Rcpp::stop("P must be square.");
  }

  const int m_targets = targets.size();
  if (m_targets == 0) {
    return Rcpp::NumericMatrix(0, 0);
  }

  for (int jj = 0; jj < m_targets; ++jj) {
    const int t = targets[jj];
    if (t < 1 || t > n) {
      Rcpp::stop("targets contains an out-of-range index.");
    }
  }

  // A0 = I - P
  SpMat A0 = SpMat(P);
  {
    double* x = A0.valuePtr();
    const int nnz = A0.nonZeros();
    for (int k = 0; k < nnz; ++k) x[k] = -x[k];
  }
  for (int i = 0; i < n; ++i) {
    A0.coeffRef(i, i) += 1.0;
  }
  A0.makeCompressed();

  // Record row -> valuePtr positions so we can overwrite one boundary row
  // without changing sparsity pattern.
  std::vector<std::vector<int>> row_slots(n);
  std::vector<int> diag_pos(n, -1);

  const int* Ap = A0.outerIndexPtr();
  const int* Ai = A0.innerIndexPtr();
  for (int col = 0; col < n; ++col) {
    for (int k = Ap[col]; k < Ap[col + 1]; ++k) {
      const int row = Ai[k];
      row_slots[row].push_back(k);
      if (row == col) diag_pos[row] = k;
    }
  }
  for (int i = 0; i < n; ++i) {
    if (diag_pos[i] < 0) {
      Rcpp::stop("Internal error: missing diagonal entry in I - P.");
    }
  }

  const std::vector<double> base_values(A0.valuePtr(), A0.valuePtr() + A0.nonZeros());

  const int out_nrow = targets_in_rows ? m_targets : n;
  const int out_ncol = targets_in_rows ? n : m_targets;
  Rcpp::NumericMatrix out(out_nrow, out_ncol);
  std::fill(out.begin(), out.end(), NA_REAL);
  double* out_ptr = out.begin();

  std::atomic<bool> has_error(false);
  std::string error_msg;

#ifdef _OPENMP
  const int use_threads = (n_threads > 0) ? n_threads : omp_get_max_threads();

  #pragma omp parallel num_threads(use_threads)
  {
    SpMat A = A0;  // thread-local numeric copy
    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(A);  // once per thread; pattern is fixed

    #pragma omp for schedule(dynamic, 1)
    for (int jj = 0; jj < m_targets; ++jj) {
      if (has_error.load()) continue;

      const int t = targets[jj] - 1;
      std::string err = compute_one_target_flighttime(
        t, jj, A, row_slots, diag_pos, base_values, solver,
        out_ptr, out_nrow, n, m_targets, targets_in_rows, hit_tol
      );

      if (!err.empty()) {
        #pragma omp critical
        {
          if (!has_error.load()) {
            has_error.store(true);
            error_msg = err;
          }
        }
      }
    }
  }
#else
  {
    SpMat A = A0;
    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(A);

    for (int jj = 0; jj < m_targets; ++jj) {
      const int t = targets[jj] - 1;

      std::string err = compute_one_target_flighttime(
        t, jj, A, row_slots, diag_pos, base_values, solver,
        out_ptr, out_nrow, n, m_targets, targets_in_rows, hit_tol
      );

      if (!err.empty()) {
        Rcpp::stop(err);
      }

      if ((jj % 32) == 0) {
        Rcpp::checkUserInterrupt();
      }
    }
  }
#endif

  if (has_error.load()) {
    Rcpp::stop(error_msg);
  }

  return out;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix cmfpt_sparse_targets_cpp(
    const Eigen::MappedSparseMatrix<double> &P,
    Rcpp::IntegerVector targets,
    bool targets_in_rows = false,
    double hit_tol = 1e-14,
    int interrupt_every = 32) {
  using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
  using Vec   = Eigen::VectorXd;

  const int n = P.rows();
  if (P.cols() != n) {
    Rcpp::stop("P must be square.");
  }

  const int m = targets.size();
  if (m == 0) {
    return Rcpp::NumericMatrix(0, 0);
  }

  // A_base = I - P
  SpMat A = P;
  {
    double *Ax0 = A.valuePtr();
    const int nnz0 = A.nonZeros();
    for (int k = 0; k < nnz0; ++k) {
      Ax0[k] = -Ax0[k];
    }
  }

  // Add identity
  for (int i = 0; i < n; ++i) {
    A.coeffRef(i, i) += 1.0;
  }
  A.makeCompressed();

  const int nnzA = A.nonZeros();
  double *Ax = A.valuePtr();
  const int *Ap = A.outerIndexPtr();
  const int *Ai = A.innerIndexPtr();

  // Cache storage positions for each row, plus diagonal position.
  // This lets us impose the boundary row x_t = const without changing the pattern.
  std::vector<std::vector<int>> row_slots(n);
  std::vector<int> diag_pos(n, -1);

  for (int col = 0; col < n; ++col) {
    for (int k = Ap[col]; k < Ap[col + 1]; ++k) {
      const int row = Ai[k];
      row_slots[row].push_back(k);
      if (row == col) {
        diag_pos[row] = k;
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    if (diag_pos[i] < 0) {
      Rcpp::stop("Missing diagonal entry after building I - P.");
    }
  }

  // Keep original values so the modified boundary row can be restored each iteration.
  std::vector<double> base_values(Ax, Ax + nnzA);

  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(A);
  
  Rcpp::NumericMatrix out = targets_in_rows
      ? Rcpp::NumericMatrix(m, n)
      : Rcpp::NumericMatrix(n, m);

  Vec rhs1(n), rhs2(n), h(n), u(n);

  for (int jj = 0; jj < m; ++jj) {
    const int t = targets[jj] - 1;
    if (t < 0 || t >= n) {
      Rcpp::stop("targets contains an out-of-range index.");
    }

    if (interrupt_every > 0 && (jj % interrupt_every == 0)) {
      Rcpp::checkUserInterrupt();
    }

    // Impose boundary row t:
    //   h_t = 1  -> row t becomes e_t^T
    //   u_t = 0  -> same matrix, different RHS
    //
    // Keep pattern unchanged by zeroing stored entries rather than removing them.
    for (int pos : row_slots[t]) {
      Ax[pos] = 0.0;
    }
    Ax[diag_pos[t]] = 1.0;

    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
      Rcpp::stop(std::string("SparseLU factorize failed at target ") + std::to_string(t + 1));
    }

    // M_t h = e_t
    rhs1.setZero();
    rhs1[t] = 1.0;
    h = solver.solve(rhs1);
    if (solver.info() != Eigen::Success) {
      Rcpp::stop(std::string("Solve for h failed at target ") + std::to_string(t + 1));
    }

    // M_t u = h, with boundary u_t = 0
    rhs2 = h;
    rhs2[t] = 0.0;
    u = solver.solve(rhs2);
    if (solver.info() != Eigen::Success) {
      Rcpp::stop(std::string("Solve for u failed at target ") + std::to_string(t + 1));
    }

    if (targets_in_rows) {
      for (int i = 0; i < n; ++i) {
        double val;
        if (i == t) {
          val = 0.0;
        } else {
          const double hi = h[i];
          val = (!std::isfinite(hi) || std::abs(hi) <= hit_tol)
              ? NA_REAL
              : (u[i] / hi);
        }
        out(jj, i) = val;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        double val;
        if (i == t) {
          val = 0.0;
        } else {
          const double hi = h[i];
          val = (!std::isfinite(hi) || std::abs(hi) <= hit_tol)
              ? NA_REAL
              : (u[i] / hi);
        }
        out(i, jj) = val;
      }
    }

    // Restore row t back to A_base before the next target.
    for (int pos : row_slots[t]) {
      Ax[pos] = base_values[pos];
    }
  }
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
// /*** R */

