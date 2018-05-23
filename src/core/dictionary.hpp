#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <iostream>
#include <string>
#include <map>

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
void dictionary::value<std::string>::print(std::ostream& stream) const;

template<>
void dictionary::value<bool>::print(std::ostream& stream) const;


#endif /* DICTIONARY_H */
