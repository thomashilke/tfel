#ifndef _LINEAR_SOLVER_H_
#define _LINEAR_SOLVER_H_

#include <numeric>
#include <vector>
#include <map>

#include <petscksp.h>

#include <spikes/array.hpp>


enum class solver {
  lapack, petsc
};

enum class method {
  dense_lu, gmres
};

enum class preconditioner {
  id, ilu
};

namespace linear_solver_impl {
  class solver_base {
  public:
    virtual ~solver_base() {}
    virtual void preallocate(const int* nz) = 0;
    virtual void add_value(std::size_t i, std::size_t j, double v) = 0;
    virtual void add_row(int i, std::size_t n, int* js, double* vs) = 0;
    virtual void assemble() = 0;
    virtual void show() = 0;
    virtual array<double> solve(const array<double>& rhs) = 0;
  };

  namespace petsc {
    struct global_initialize {
      static global_initialize& instance() {
	if (not inst)
	  inst = new global_initialize();
	return *inst;
      }

      static void release() {
	delete inst;
	inst = nullptr;
      }
      
    private:
      global_initialize() {
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
      }
      ~global_initialize() {
	PetscFinalize();
      }

    private:
      static global_initialize* inst;
    };
  }
  
  class petsc_gmres_ilu: public solver_base {
  public:
    petsc_gmres_ilu(std::size_t problem_size): verbose(false) {
      PetscErrorCode ierr;

      petsc::global_initialize::instance();
      
      ierr = MatCreate(PETSC_COMM_WORLD, &A);CHKERRV(ierr);
      ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, problem_size, problem_size);CHKERRV(ierr);
      ierr = MatSetType(A, MATSEQAIJ);CHKERRV(ierr);

      ierr = VecCreate(PETSC_COMM_WORLD, &b);CHKERRV(ierr);
      ierr = VecSetSizes(b, PETSC_DECIDE, problem_size);CHKERRV(ierr);
      ierr = VecSetFromOptions(b);CHKERRV(ierr);

      ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRV(ierr);
      ierr = KSPSetOperators(ksp, A, A);CHKERRV(ierr);
      ierr = KSPSetType(ksp,KSPGMRES);CHKERRV(ierr);
      //ierr = KSPSetType(ksp,KSPPREONLY);CHKERRV(ierr);
      //ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);CHKERRV(ierr);

      ierr = KSPGetPC(ksp, &pc);CHKERRV(ierr);

      //ierr = PCSetType(pc,PCLU);CHKERRV(ierr);
      ierr = PCSetType(pc,PCILU);CHKERRV(ierr);
      ierr = PCFactorSetLevels(pc, 2);CHKERRV(ierr);

      //ierr = PCFactorSetAllowDiagonalFill(pc, PETSC_TRUE);CHKERRV(ierr);
      //ierr = PCFactorSetMatOrderingType(pc, MATORDERINGRCM);CHKERRV(ierr);
      

      ierr = KSPSetTolerances(ksp, 1.e-8, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRV(ierr);
      //ierr = KSPGMRESSetOrthogonalization(ksp, KSPGMRESClassicalGramSchmidtOrthogonalization);CHKERRV(ierr);

      ierr = PetscViewerCreate(PETSC_COMM_WORLD, &v);CHKERRV(ierr);
      ierr = PetscViewerSetType(v, PETSCVIEWERASCII);CHKERRV(ierr);
      ierr = PetscViewerPushFormat(v, PETSC_VIEWER_ASCII_MATLAB);CHKERRV(ierr);
      ierr = PetscViewerAndFormatCreate(v, PETSC_VIEWER_DEFAULT, &vaf);CHKERRV(ierr);
      if (verbose)
	ierr = KSPMonitorSet(ksp,
			     reinterpret_cast<PetscErrorCode (*)(KSP,
								 PetscInt,
								 PetscReal,
								 void *)>(KSPMonitorTrueResidualNorm),
			     vaf,
			     nullptr);CHKERRV(ierr);
      ierr = KSPSetFromOptions(ksp);CHKERRV(ierr);
    }

    ~petsc_gmres_ilu() {
      PetscErrorCode ierr;
      
      ierr = MatDestroy(&A);CHKERRV(ierr);
      ierr = VecDestroy(&x);CHKERRV(ierr);
      ierr = VecDestroy(&b);CHKERRV(ierr);
      ierr = KSPDestroy(&ksp);CHKERRV(ierr);
    }

    virtual void preallocate(const int* nz) {
      PetscErrorCode ierr;

      ierr = MatSeqAIJSetPreallocation(A, 0, nz);CHKERRV(ierr);
      ierr = MatSetUp(A);CHKERRV(ierr);
    }

    void add_value(std::size_t i, std::size_t j, double v) {
      PetscErrorCode ierr;
      ierr = MatSetValue(A, i, j, v, ADD_VALUES);CHKERRV(ierr);
    }

    void add_row(int i, std::size_t n, int* js, double* vs) {
      PetscErrorCode ierr;
      ierr = MatSetValues(A, 1, &i, n, js, vs, ADD_VALUES);CHKERRV(ierr);
    }
    
    void assemble() {
      PetscErrorCode ierr;
      ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
      ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
    }

    void show() {
      PetscErrorCode ierr;
      ierr = MatView(A, v);CHKERRV(ierr);
    }
    
    virtual array<double> solve(const array<double>& rhs) {
      PetscErrorCode ierr;
      for (std::size_t i(0); i < rhs.get_size(0); ++i) {
	ierr = VecSetValue(b, i, rhs.at(i), INSERT_VALUES);CHKERRCONTINUE(ierr);
      }
      ierr = VecAssemblyBegin(b);CHKERRCONTINUE(ierr);
      ierr = VecAssemblyBegin(b);CHKERRCONTINUE(ierr);

      ierr = VecDuplicate(b, &x);CHKERRCONTINUE(ierr);
      ierr = KSPSolve(ksp, b, x);CHKERRCONTINUE(ierr);

      array<double> result{rhs.get_size(0)};
      std::vector<PetscInt> ix(rhs.get_size(0));
      std::iota(ix.begin(), ix.end(), 0);

      ierr = VecGetValues(x, ix.size(), &ix[0], &result.at(0));CHKERRCONTINUE(ierr);
      return result;
    }
    
  private:
    const bool verbose;
    
    Mat A;
    Vec x, b;

    KSP ksp;
    PC pc;

    PetscViewer v;
    PetscViewerAndFormat* vaf;
  };

  class lapack_dense_lu: public solver_base {
  public:
    lapack_dense_lu(std::size_t problem_size)
      : a{problem_size, problem_size} {
      double* data(a.get_data());
      std::fill(data, data + problem_size * problem_size, 0.0);
    }

    virtual void preallocate(const int* nz) {}
	  
    void add_value(std::size_t i, std::size_t j, double v) {
      a.at(i, j) = v;
    }

    void add_row(int i, std::size_t n, int* js, double* vs) {}

    void assemble() {}
    void show() {}
    virtual array<double> solve(const array<double>& /*rhs*/) {
      throw std::string("unimplemented");
    }

  private:
    array<double> a;
  };
};

class linear_solver {
public:
  linear_solver() {
    solver_factories[solver::lapack] = build_solver<linear_solver_impl::lapack_dense_lu>;
    solver_factories[solver::petsc] = build_solver<linear_solver_impl::petsc_gmres_ilu>;
  }
  
  linear_solver_impl::solver_base* get_solver(solver s, method m, preconditioner pc, std::size_t problem_size) {
    return solver_factories[s](m, pc, problem_size);
  }

private:
  typedef linear_solver_impl::solver_base* (*solver_factory)(method, preconditioner, std::size_t);
  std::map<solver, solver_factory> solver_factories;
  
  template<typename S>
  static linear_solver_impl::solver_base* build_solver(method m, preconditioner pc, std::size_t problem_size) {
    return new S(problem_size);
  }
};


#endif /* _LINEAR_SOLVER_H_ */
