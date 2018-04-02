
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <vector>

#include <spikes/array.hpp>

#include <lapacke.h>
#include <petscksp.h>
#include <petsc/private/kspimpl.h>
#include <petscdmshell.h>


class dictionary {
private:
  struct basic_value {
    virtual void print(std::ostream& stream) const = 0;
  };

  template<typename T>
  struct value: public basic_value {
    value(const T& v): v(v) {}
    
    virtual void print(std::ostream& stream) const;
    
    T v;
  };

  using value_ptr = std::shared_ptr<basic_value>;
  
public:
  dictionary() {}

  template<typename T>
  dictionary& set(const std::string& key, const T& v) {
    auto p = new value<T>(v);
    kv[key] = value_ptr(p);
    return *this;
  }

  template<std::size_t n>
  dictionary& set(const std::string& key, const char (&str)[n]) {
    auto p = new value<std::string>(str);
    kv[key] = value_ptr(p);
    return *this;
  }

  template<typename T>
  const T& get(const std::string& key) const {
    auto item(kv.find(key));
    if (item == kv.end())
      throw std::string("dictionary::get: key '") + key + "' not found.";

    const value<T>* p(dynamic_cast<const value<T>*>((item->second).get()));
    if (p == nullptr)
      throw std::string("dictionary::get: key '") + key + "' is of wrong type.";

    return p->v;
  }

  bool key_exists(const std::string& key) const {
    return kv.find(key) != kv.end();
  }

  template<typename Iterator>
  bool keys_exist(Iterator begin, Iterator end) const {
    while (begin != end) {
      if (not key_exists(*begin))
        return false;
      ++begin;
    }

    return true;
  }

  void clear() {
    kv.clear();
  }

  void print(std::ostream& stream) const {
    bool first_item(true);

    stream << "dictionary {";
    for (const auto& item: kv) {
      if (not first_item)
        stream << "," << std::endl;
      else
        stream << std::endl;

      stream << "  \"" << item.first << "\": ";
      item.second->print(stream);
      
      first_item = false;
    }
    stream << std::endl << "}" << std::endl;
  }

private:  
  std::map<std::string, value_ptr> kv;
};

template<typename T>
void dictionary::value<T>::print(std::ostream& stream) const {
  stream << v;
}

template<>
void dictionary::value<std::string>::print(std::ostream& stream) const {
  stream << "\"" << v << "\"";
}

template<>
void dictionary::value<bool>::print(std::ostream& stream) const {
  stream << std::boolalpha << v;
}

class matrix;
class sparse_matrix;
class dense_matrix;
class crs_matrix;
class skyline_matrix;

namespace solver {
  class basic_solver {
  public:
    virtual ~basic_solver() {}

    virtual void set_operator(const matrix& m) = 0;
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
      
      void set_operator(const matrix& m) {
        
      }
      
      void set_operator(const sparse_matrix& m) {

      }
      
      void set_operator(const dense_matrix& m) {

      }

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
          "rtol",   "abstol",
          "dtol",   "ilufill"
        };
  
        if (not params.keys_exist(expected_keys.begin(), expected_keys.end()))
          throw std::string("solver::petsc::gmres_ilu: missing key(s) "
                            "in parameter dictionary.");
  
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

      ~gmres_ilu() {
        PetscErrorCode ierr;
        
        ierr = MatDestroy(&a);CHKERRV(ierr);
        ierr = VecDestroy(&b);CHKERRV(ierr);
        ierr = KSPDestroy(&ksp);CHKERRV(ierr);
      }
      
      void set_operator(const matrix& m);

      void set_operator(const sparse_matrix& m) {

      }
      
      void set_operator(const dense_matrix& m) {

      }
      
      bool solve(const array<double>& rhs,
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
                  << " resid norm " << rnorm
                  << " true resid norm " << truenorm
                  << "||r(i)||/||b||" << truenorm/bnorm
                  << std::endl;
        
        return 0;
      }
    };
  }
}

solver::petsc::initialize* solver::petsc::initialize::inst(nullptr);

class matrix {
public:
  virtual ~matrix() {}

  virtual void populate_solver(solver::basic_solver* s) const = 0;
  
  virtual std::size_t get_row_number() const = 0;
  virtual std::size_t get_column_number() const = 0;

  virtual std::size_t get_nz_element_number() const = 0;

  virtual void fill(double v) = 0;
  
  virtual void set(std::size_t i, std::size_t j, double v) = 0;
  virtual void add(std::size_t i, std::size_t j, double v) = 0;
  virtual double get(std::size_t i, std::size_t j) const = 0;
};

class sparse_matrix: public matrix {
public:
  sparse_matrix(std::size_t n_row, std::size_t n_column)
    : n_row(n_row), n_column(n_column) {}

  virtual void populate_solver(solver::basic_solver* s) const {
    s->set_operator(*this);
  }
  
  virtual std::size_t get_row_number() const { return n_row; }
  virtual std::size_t get_column_number() const { return n_column; }

  virtual std::size_t get_nz_element_number() const {
    return values.size();
  }
  
  virtual void fill(double v) {
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
  
private:
  std::size_t n_row, n_column;
  std::map<std::pair<std::size_t, std::size_t>, double> values;
};

class dense_matrix: public matrix {
public:
  dense_matrix(std::size_t n_row, std::size_t n_column): values{n_row, n_column} {}

  virtual void populate_solver(solver::basic_solver* s) const {
    s->set_operator(*this);
  }
  
  virtual std::size_t get_row_number() const { return values.get_size(0); }
  virtual std::size_t get_column_number() const { return values.get_size(1); }

  virtual std::size_t get_nz_element_number() const {
    return values.get_size(0) * values.get_size(1);
  }
  
  virtual void fill(double v) { values.fill(v); }
  
  virtual void set(std::size_t i, std::size_t j, double v) {
    values.at(i, j) = v;
  }
    
  virtual void add(std::size_t i, std::size_t j, double v)  {
    values.at(i, j) += v;
  }
  
  virtual double get(std::size_t i, std::size_t j) const {
    return values.at(i, j);
  }
  
private:
  array<double> values;
};


bool solver::petsc::gmres_ilu::solve(const array<double>& rhs,
                                     array<double>& x,
                                     dictionary& report) {
  PetscErrorCode ierr;
  
  ierr = VecZeroEntries(b);
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

/*
void solver::petsc::gmres_ilu::set_operator(const matrix& m) {
  PetscErrorCode ierr;
  
  ierr = MatZeroEntries(a);CHKERRV(ierr);
  
  ierr = MatSetSizes(a,
                     m.get_row_number(), m.get_column_number(),
                     m.get_row_number(), m.get_column_number());CHKERRV(ierr);
  ierr = VecSetSizes(b,
                     m.get_row_number(), m.get_column_number());CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(a, m.get_nz_element_number(), nullptr);CHKERRV(ierr);
  ierr = VecSetUp(b);CHKERRV(ierr);
  
  m.populate_solver(this);

  ierr = MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
  ierr = KSPSetOperators(ksp, a, a);CHKERRV(ierr);
}

void solver::lapack::lu::set_operator(const matrix& m) {
  report.clear();
  
  m.populate_solver(this);
  pivots = array<lapack_int>{m.get_row_number()};

  lapack_int n(m.get_row_number());
  lapack_int info(LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
                                 n, n, data.get_data(), n,
                                 pivots.get_data()));

  if (info > 0) {
    this->report.set("error", "singular matrix");
  } else if (info < 0) {
    this->report.set("error", "dgetrf invalid parameter");
  }

  valid_decomposition = true;
}
*/

void fd_heat(sparse_matrix& m, array<double>& rhs, std::size_t n) {
  for (std::size_t i(0); i < n; ++i) {
    m.set(i, i, 1.0);
    rhs.at(i) = 1.0;
  }
}

int main() {
  try {
    sparse_matrix a(10, 10);
    array<double> rhs{10};
    array<double> x{10};

    fd_heat(a, rhs, 10);
    
    dictionary p(dictionary()
                 .set("maxits",  2000u)
                 .set("restart", 1000u)
                 .set("rtol",    1.e-8)
                 .set("abstol",  1.e-50)
                 .set("dtol",    1.e20)
                 .set("ilufill", 2u)
                 .set("test_bool", true)
                 .set("test_string", "str"));
    p.print(std::cout);
    
    //solver::petsc::gmres_ilu s(p);
    solver::lapack::lu s;
    

    s.set_operator(a);
    dictionary report;
    if (s.solve(rhs, x, report)) {
      std::cout << "success." << std::endl;
      //std::cout << report.get<unsigned int>("iterations") << " iterations." << std::endl;
    } else {
      report.print(std::cout);
      std::cout << "failed to converge: " << report.get<std::string>("error") << std::endl;
    }
  
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  
  return 0;
}
