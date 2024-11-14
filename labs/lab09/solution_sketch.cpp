class HeatEquation
{
public:
  HeatEquation();

  void
  run();

private:
  void
  setup_system();

  void
  assemble_system();

  double
  get_source(const double time, const Point<2> &point) const;

  Vector<double>
  evaluate_diffusion(const double time, const Vector<double> &y) const;

  void
  output_results(const double                     time,
                 const unsigned int               time_step,
                 TimeStepping::runge_kutta_method method) const;

  void
  explicit_method(const TimeStepping::runge_kutta_method method,
                  const unsigned int                     n_time_steps,
                  const double                           initial_time,
                  const double                           final_time);

  const unsigned int fe_degree;

  const double diffusion_coefficient;

  Triangulation<2>          triangulation;
  const FE_Q<2>             fe;
  DoFHandler<2>             dof_handler;
  AffineConstraints<double> constraint_matrix;
  SparsityPattern           sparsity_pattern;

  SparseMatrix<double> system_matrix;
  SparseMatrix<double> mass_matrix;
  SparseDirectUMFPACK  inverse_mass_matrix;

  Vector<double> solution;
};

void
HeatEquation::assemble_system()
{
  // Assemble system matrix (=-S)
  // Assemble mass matrix (=M)

  // Initialize solver to invert M (inverse_mass_matrix.initialize(mass_matrix).
}

Vector<double>
HeatEquation::evaluate_diffusion(const double time, const Vector<double> &y)
{
  Vector<double> tmp(n_dofs);
  tmp = 0;
  system_matrix.vmult(tmp, y);

  // Compute F(time) and add contribution to tmp.
  Vector<double> value;
  inverse_mass_matrix.vmult(value, tmp);

  return value;
}

void
HeatEquation::run()
{
  const double time_step =
    (final_time - initial_time) / static_cast<double>(n_time_steps);

  solution = 0.;
  constraint_matrix.distribute(solution);

  TimeStepping::ExplicitRungeKutta<Vector<double>> explicit_runge_kutta(TimeStepping::runge_kutta_method::RK_CLASSIC_FOURTH_ORDER);
  output_results(initial_time, 0, method);

  DiscreteTime time(initial_time, final_time, time_step);

  while (time.is_at_end() == false)
    {
      explicit_runge_kutta.evolve_one_time_step(
        [this](const double time, const Vector<double> &y) {
          return this->evaluate_diffusion(time, y);
        },
        time.get_current_time(),
        time.get_next_step_size(),
        solution);
      time.advance_time();

      constraint_matrix.distribute(solution);

      if (time.get_step_number() % 10 == 0)
        output_results(time.get_current_time(), time.get_step_number(), method);
    }
}
