#ifndef IMPORTER_H
#define IMPORTER_H

#include <alucelldb/alucelldb.hpp>


#include "../core/cell.hpp"
#include "../core/mesh.hpp"
#include "../core/mesh_data.hpp"


#ifdef ENABLE_ALUCELL
#include "alucell_importer.hpp"
#endif

#ifdef ENABLE_FREEFEM
#include "freefem_importer.hpp"
#endif

#ifdef ENABLE_GMSH
#include "gmsh_importer.hpp"
#endif

#endif /* IMPORTER_H */
