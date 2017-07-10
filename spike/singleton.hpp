#ifndef SINGLETON_H
#define SINGLETON_H

template<typename T>
class singleton {
public:
  static T& instance() {
    if (not inst)
      inst = new T();
    return inst;
  }
  static void release() {
    delete inst;
    inst = nullptr;
  }
protected:
  virtual ~singleton() {}
  
private:
  static T* inst;
};

template<typename T>
T* singleton<T>::inst(nullptr);

#endif /* SINGLETON_H */
