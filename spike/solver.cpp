
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <vector>

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

namespace solver {
  namespace petsc {
    int gmres_ilu(const dictionary& params) {
      std::vector<std::string> expected_keys {
        "maxits", "restart",
        "rtol",   "abstol",
        "dtol",   "ilufill"
      };

      if (not params.keys_exist(expected_keys.begin(), expected_keys.end()))
        throw std::string("solver::petsc::gmres_ilu: missing key(s) "
                          "in parameter dictionary.");

      params.get<unsigned int>("maxits");
      params.get<unsigned int>("restart");
      params.get<double>("rtol");
      params.get<double>("abstol");
      params.get<double>("dtol");
      params.get<unsigned int>("ilufill");
      params.get<bool>("test_bool");
      params.get<std::string>("test_string");
      
      return 1;
    }
  }
}

int main() {
  try {
  auto stat_stokes_solver_params(dictionary()
                                 .set("maxits",  2000u)
                                 .set("restart", 1000u)
                                 .set("rtol",    1.e-8)
                                 .set("abstol",  1.e-50)
                                 .set("dtol",    1.e20)
                                 .set("ilufill", 2u)
                                 .set("test_bool", true)
                                 .set("test_string", "str"));
  
  auto stat_stokes_solver = solver::petsc::gmres_ilu(stat_stokes_solver_params);

  stat_stokes_solver_params.print(std::cout);

  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  
  return 0;
}
