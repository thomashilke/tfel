#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>

#include <spikes/array.hpp>

#include <lapacke.h>
#include <petscksp.h>
#include <petsc/private/kspimpl.h>
#include <petscdmshell.h>

#include "meta.hpp"
#include "dictionary.hpp"


class matrix;
class sparse_matrix;
class dense_matrix;
class crs_matrix;
class skyline_matrix;

namespace solver {
  class basic_solver {
  public:
    virtual ~basic_solver() {}

    void set_operator(const matrix& m);
    virtual void set_operator(const sparse_matrix& m) = 0;
    virtual void set_operator(const dense_matrix& m) = 0;
    
    virtual bool solve(const array<double>& rhs,
                       array<double>& x,
                       dictionary& report) = 0;
  };

  namespace lapack {
    class lu: public basic_solver {
    public:
      lu(): data{1}, pivots{1}, valid_decomposition(false) {}
      virtual ~lu() {}
      
      void set_operator(const sparse_matrix& m);
      void set_operator(const dense_matrix& m);

      void set_operator_size(std::size_t n) {
        data = array<double>{n, n};
      }

      void set(std::size_t i, std::size_t j, double v) {
        data.at(i, j) = v;
      }

      bool solve(const array<double>& rhs,
                 array<double>& x,
                 dictionary& report) {
        if (not valid_decomposition) {
          report = this->report;
          return false;
        }

        x = rhs;
        
        lapack_int n(data.get_size(0));
        lapack_int info(LAPACKE_dgetrs(LAPACK_ROW_MAJOR,
                                       'N',
                                       n, 1, data.get_data(), n, pivots.get_data(),
                                       x.get_data(), n));

        if (info < 0) {
          report.set("error", "dgetrs invalid parameter");
          return false;
        }

        return true;
      }

    private:
      array<double> data;
      array<lapack_int> pivots;
      
      bool valid_decomposition;
      dictionary report;

    private:
      void do_lu_decomposition();
    };
  }
  
  namespace petsc {
    class initialize {
    public:
      static initialize& instance() {
        if (not inst)
          inst = new initialize();
        return *inst;
      }
      
      static void release() {
        delete inst;
        inst = nullptr;
      }

    private:
      initialize() {
        PetscInitialize(nullptr, nullptr, nullptr, nullptr);
      }
      
      ~initialize() {
        PetscFinalize();
      }

      static initialize* inst;
    };
    
    
    class gmres_ilu: public basic_solver {
    public:
      gmres_ilu(const dictionary& params) {
        std::vector<std::string> expected_keys {
          "maxits", "restart",
          "rtol",   "atol",
          "dtol",   "ilufill",
        };
  
        if (not params.keys_exist(expected_keys.begin(), expected_keys.end()))
          throw std::string("solver::petsc::gmres_ilu: missing key(s) "
                            "in parameter dictionary.");
        
        const auto& petsc_init(petsc::initialize::instance());
        ignore_unused(petsc_init);
        
        PetscErrorCode ierr;
        ierr = MatCreate(PETSC_COMM_WORLD, &a);CHKERRV(ierr);
        ierr = MatSetType(a, MATSEQAIJ);CHKERRV(ierr);

        ierr = VecCreate(PETSC_COMM_WORLD, &b);CHKERRV(ierr);
        ierr = VecSetType(b, VECSEQ);CHKERRV(ierr);

        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRV(ierr);

        
        ierr = KSPSetType(ksp, KSPGMRES);CHKERRV(ierr);
        ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);CHKERRV(ierr);
        ierr = KSPSetTolerances(ksp,
                                params.get<double>("rtol"),
                                params.get<double>("atol"),
                                params.get<double>("dtol"),
                                params.get<unsigned int>("maxits"));CHKERRV(ierr);
        ierr = KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);CHKERRV(ierr);
        ierr = KSPGMRESSetRestart(ksp, params.get<unsigned int>("restart"));CHKERRV(ierr);

        ierr = KSPMonitorSet(ksp, monitor, nullptr, nullptr);CHKERRV(ierr);

        PC pc;
        ierr = KSPGetPC(ksp, &pc);CHKERRV(ierr);
        ierr = PCSetType(pc, PCILU);CHKERRV(ierr);
        ierr = PCFactorSetLevels(pc, params.get<unsigned int>("ilufill"));CHKERRV(ierr); 
      }

      virtual ~gmres_ilu() {
        PetscErrorCode ierr;
        
        ierr = MatDestroy(&a);CHKERRV(ierr);
        ierr = VecDestroy(&b);CHKERRV(ierr);
        ierr = KSPDestroy(&ksp);CHKERRV(ierr);
      }
      
      virtual void set_operator(const sparse_matrix& m);
      virtual void set_operator(const dense_matrix& m);
      
      virtual bool solve(const array<double>& rhs,
                         array<double>& x,
                         dictionary& report);
      
    private:
      Mat a;
      Vec b;
      KSP ksp;

      static PetscErrorCode monitor(KSP ksp, PetscInt it, PetscReal rnorm, void*) {
        Vec            resid;
        PetscReal      truenorm,bnorm;
        char           normtype[256];
        
        if (it == 0 && ((PetscObject)ksp)->prefix)
          std::cout << "Residual norms for " << ((PetscObject)ksp)->prefix << " solve." << std::endl;

        KSPBuildResidual(ksp,NULL,NULL,&resid);
        VecNorm(resid,NORM_2,&truenorm);
        VecDestroy(&resid);
        VecNorm(ksp->vec_rhs,NORM_2,&bnorm);
        
        PetscStrncpy(normtype,KSPNormTypes[ksp->normtype],sizeof(normtype));
        PetscStrtolower(normtype);
        
        std::cout << it << " KSP " << normtype
                  << " resid norm " << std::setw(14) << std::right << rnorm
                  << " true resid norm " << std::setw(14) << std::right << truenorm
                  << "||r(i)||/||b||" << std::setw(14) << std::right << truenorm/bnorm
                  << std::endl;
        
        return 0;
      }
    };
  }
}



class matrix {
public:
  virtual ~matrix() {}

  virtual void populate_solver(solver::basic_solver& s) const = 0;
  
  virtual std::size_t get_row_number() const = 0;
  virtual std::size_t get_column_number() const = 0;

  virtual std::size_t get_nz_element_number() const = 0;

  virtual void clear() = 0;
  
  virtual void set(std::size_t i, std::size_t j, double v) = 0;
  virtual void add(std::size_t i, std::size_t j, double v) = 0;
  virtual double get(std::size_t i, std::size_t j) const = 0;
  virtual double& get(std::size_t i, std::size_t j) = 0;
};


class sparse_matrix: public matrix {
public:
  friend class solver::lapack::lu;
  friend class solver::petsc::gmres_ilu;
  
  sparse_matrix(std::size_t n_row, std::size_t n_column)
    : n_row(n_row), n_column(n_column) {}

  virtual void populate_solver(solver::basic_solver& s) const {
    s.set_operator(*this);
  }
  
  virtual std::size_t get_row_number() const { return n_row; }
  virtual std::size_t get_column_number() const { return n_column; }

  virtual std::size_t get_nz_element_number() const {
    return values.size();
  }
  
  virtual void clear() {
    values.clear();
  }
  
  virtual void set(std::size_t i, std::size_t j, double v) {
    values[std::make_pair(i, j)] = v;
  }
  
  virtual void add(std::size_t i, std::size_t j, double v) {
    values[std::make_pair(i, j)] += v;
  }
  
  virtual double get(std::size_t i, std::size_t j) const {
    auto item(values.find(std::make_pair(i,j)));
    if (item == values.end())
      return 0.0;
    else
      return item->second;
  }

  virtual double& get(std::size_t i, std::size_t j) {
    auto result(values.insert(std::make_pair(std::make_pair(i, j), 0.0)));
    return result.first->second;
  }
  
private:
  std::size_t n_row, n_column;
  std::map<std::pair<std::size_t, std::size_t>, double> values;
};


class dense_matrix: public matrix {
public:
  friend class solver::lapack::lu;
  friend class solver::petsc::gmres_ilu;
  
  dense_matrix(std::size_t n_row, std::size_t n_column): values{n_row, n_column} {}

  virtual void populate_solver(solver::basic_solver& s) const {
    s.set_operator(*this);
  }
  
  virtual std::size_t get_row_number() const { return values.get_size(0); }
  virtual std::size_t get_column_number() const { return values.get_size(1); }

  virtual std::size_t get_nz_element_number() const {
    return values.get_size(0) * values.get_size(1);
  }
  
  virtual void clear() { values.fill(0.0); }
  
  virtual void set(std::size_t i, std::size_t j, double v) {
    values.at(i, j) = v;
  }
    
  virtual void add(std::size_t i, std::size_t j, double v)  {
    values.at(i, j) += v;
  }
  
  virtual double get(std::size_t i, std::size_t j) const {
    return values.at(i, j);
  }

  virtual double& get(std::size_t i, std::size_t j) {
    return values.at(i, j);
  }
  
private:
  array<double> values;
};


#endif /* SOLVER_H */
