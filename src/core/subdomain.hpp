#ifndef SUBDOMAIN_H
#define SUBDOMAIN_H

namespace cell {
  namespace detail {

    template<std::size_t n>
    struct subdomain {
      subdomain(): fill(0) {}

      void insert(unsigned int node) {
	if (fill and nodes[fill - 1] >= node)
	  throw std::string("nodes must be ordered when inserted into a subdomain");
	
	nodes[fill] = node;
	fill += 1;
      }

      std::size_t size() const { return fill; }

      unsigned int* begin() { return nodes; }
      unsigned int* end() { return nodes + fill; }

      const unsigned int* begin() const { return nodes; }
      const unsigned int* end() const { return nodes + fill; }

      bool operator==(const subdomain<n>& s) const {
	if (fill != s.fill)
	  return false;
	
	for (unsigned int i(0); i < fill; ++i)
	  if (nodes[i] != s.nodes[i])
	    return false;
      }

      bool operator<(const subdomain<n>& s) const {
	if (fill != s.fill)
	  return fill < s.fill;
	for (unsigned int i(0); i < fill; ++i)
	  if (nodes[i] != s.nodes[i])
	    return nodes[i] < s.nodes[i];
	return false;
      }
      
      unsigned int nodes[n];
      unsigned short fill;
    };
  }
  
  typedef detail::subdomain<4> subdomain_type;

}

#endif /* SUBDOMAIN_H */
