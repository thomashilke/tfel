#include "solver.hpp"

void solver::basic_solver::set_operator(const matrix& m) {
  m.populate_solver(*this);
}


void solver::lapack::lu::set_operator(const sparse_matrix& m) {
  report.clear();

  /*
   *  Fill the data array
   */
  data.fill(0.0);
  for (auto v: m.values)
    data.at(v.first.first, v.first.second) = v.second;

  do_lu_decomposition();
}

void solver::lapack::lu::set_operator(const dense_matrix& m) {
  report.clear();
  
  /*
   *  Fill the data array
   */
  data = m.values;


  do_lu_decomposition();
}

void solver::lapack::lu::do_lu_decomposition() {
  /*
   *  LU decomposition
   */
  pivots = array<lapack_int>{data.get_size(0)};

  lapack_int n(data.get_size(0));
  lapack_int info(LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
                                 n, n, data.get_data(), n,
                                 pivots.get_data()));

  if (info > 0) {
    this->report.set("error", "singular matrix");
  } else if (info < 0) {
    this->report.set("error", "dgetrf invalid parameter");
  }

  valid_decomposition = info == 0;
}


solver::petsc::initialize* solver::petsc::initialize::inst = nullptr;


void solver::petsc::gmres_ilu::set_operator(const sparse_matrix& m) {
  PetscErrorCode ierr;

  ierr = MatSetSizes(a,
                     m.get_row_number(), m.get_column_number(),
                     m.get_row_number(), m.get_column_number());CHKERRV(ierr);
  
  ierr = VecSetSizes(b,
                     m.get_column_number(), m.get_column_number());CHKERRV(ierr);

  
  /*
   *  Convert sparse matrix to CRS representation
   */
  std::vector<int>
    row(m.get_row_number() + 1),
    col(m.get_nz_element_number());
  std::vector<double>
    val(m.get_nz_element_number());
    
  std::size_t row_id(0), val_id(0);
  row[0] = row_id;
  for (const auto& v: m.values) {
    while (row_id < v.first.first) {
      ++row_id;
      row[row_id] = val_id;
    }
    col[val_id] = v.first.second;
    val[val_id] = v.second;
    ++val_id;
  }
  row.back() = val_id;

  
  /*
   *  Count the number of non-zero per row
   */
  std::vector<int> nnz(m.get_row_number());
  for (std::size_t n(0); n < m.get_row_number(); ++n)
    nnz[n] = row[n + 1] - row[n];

  
  /*
   *  Set preallocation of the sequential matrix storage
   */
  ierr = MatSeqAIJSetPreallocation(a, 0, &nnz[0]);

  
  /*
   * Assemble the matrix row by row
   */
  for (int row_id(0); row_id < row.size() - 1; ++row_id) {
    if (row[row_id + 1] > row[row_id]) {
      ierr = MatSetValues(a, 1, &row_id,
                          nnz[row_id], &col[row[row_id]],
                          &val[row[row_id]],
                          INSERT_VALUES);CHKERRV(ierr);
    }
  }

  
  /*
   *  Finilize assembly
   */
  ierr = MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = KSPSetOperators(ksp, a, a);CHKERRV(ierr);  
}


void solver::petsc::gmres_ilu::set_operator(const dense_matrix& m) {
  PetscErrorCode ierr;

  ierr = MatSetSizes(a,
                     m.get_row_number(), m.get_column_number(),
                     m.get_row_number(), m.get_column_number());CHKERRV(ierr);
  
  ierr = VecSetSizes(b,
                     m.get_column_number(), m.get_column_number());CHKERRV(ierr);

  ierr = MatSeqAIJSetPreallocation(a, m.get_row_number() * m.get_column_number(),
                                   nullptr);CHKERRV(ierr);

  std::vector<int> is(m.get_column_number());
  std::iota(is.begin(), is.end(), 0);
  for (int row_id(0); row_id < m.get_row_number(); ++row_id) {
    ierr = MatSetValues(a, 1, &row_id,
                        is.size(), &is[0],
                        &m.values.at(row_id, 0),
                        INSERT_VALUES);CHKERRV(ierr);
  }

  ierr = MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = KSPSetOperators(ksp, a, a);CHKERRV(ierr);  
}


bool solver::petsc::gmres_ilu::solve(const array<double>& rhs,
                                     array<double>& x,
                                     dictionary& report) {
  PetscErrorCode ierr;
  
  ierr = VecZeroEntries(b);CHKERRCONTINUE(ierr);
  for (std::size_t i(0); i < rhs.get_size(0); ++i) {
    ierr = VecSetValue(b, i, rhs.at(i), INSERT_VALUES);CHKERRCONTINUE(ierr);
  }
  ierr = VecAssemblyBegin(b);CHKERRCONTINUE(ierr);
  ierr = VecAssemblyEnd(b);CHKERRCONTINUE(ierr);

  Vec y;
  ierr = VecDuplicate(b, &y);CHKERRCONTINUE(ierr);
  ierr = KSPSolve(ksp, b, y);CHKERRCONTINUE(ierr);

  std::vector<PetscInt> iy(rhs.get_size(0));
  std::iota(iy.begin(), iy.end(), 0);

  ierr = VecGetValues(y, iy.size(), &iy[0], &x.at(0));CHKERRCONTINUE(ierr);
  ierr = VecDestroy(&y);CHKERRCONTINUE(ierr);
  return true;
}

