#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <iostream>

using namespace dealii;

int
main()
{
  // Triangulation.
  Triangulation<2> tria;

  Point<2> center(0, 0);

  double inner_radius = 0.5, outer_radius = 1;

  GridGenerator::hyper_shell(tria, center, inner_radius, outer_radius, 5);

  // Refinement.
  for (unsigned step = 0; step < 3; ++step)
    {
      for (auto &cell : tria.active_cell_iterators())
        {
          for (auto v : cell->vertex_indices())
            {
              Point<2> vertex = cell->vertex(v);

              const double distance_from_center = vertex.distance(center);

              if (std::fabs(distance_from_center - inner_radius) <
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }
      tria.execute_coarsening_and_refinement();
    }

  // Output mesh.
  GridOut grid_out;
  grid_out.write_msh(tria, "first_mesh.msh");

  // Defining FE space.
  FE_Q<2> fe(1); // Degree 1 in 2D.

  // Distributing degrees of freedom.
  DoFHandler<2> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // Renumbering.
  DoFRenumbering::Cuthill_McKee(dof_handler);

  // Creating dynamic sparsity pattern (inefficient but memory conservative).
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);

  // Copy to static sparsity pattern (both efficient and memory conservative).
  SparsityPattern sp;
  sp.copy_from(dsp);

  // Output sparsity pattern.
  std::ofstream out("sparsity-pattern.svg");
  sp.print_svg(out);

  // Initialize system matrix and a solution vector.
  SparseMatrix<double> system_matrix;
  system_matrix.reinit(sp);

  Vector<double> vector;
  vector.reinit(dof_handler.n_dofs());

  return 0;
}
