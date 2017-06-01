#ifndef _MESH_H_
#define _MESH_H_

template<typename fe>
class mesh {
public:
  typedef fe fe_type;
  typedef typename fe::cell_type cell_type;

private:
  array<double> nodes;
  array<unsigned int> elements;
};

#endif /* _MESH_H_ */
