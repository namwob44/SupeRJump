// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppEigen.h>
#include <RcppParallel.h>

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
 
  return Rcpp::wrap(temp_output.colwise().sum());
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
// /*** R */

