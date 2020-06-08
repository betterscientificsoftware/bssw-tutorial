process_args(int argc, char **argv)
// Utilities
double
l2_norm(int n, double const *a, double const *b)

void
copy(int n, double *dst, double const *src)


void
write_array(int t, int n, double dx, double const *a)

void
set_initial_condition(int n, double *a, double dx, char const *ic)

static void
r83_np_fa(int n, double *a)

void
initialize_crankn(int n,
    double alpha, double dx, double dt,
    double **_cn_Amat)
static void 
r83_np_sl ( int n, double const *a_lu, double const *b, double *x)

bool
update_solution_crankn(int n,
    double *curr, double const *last,
    double const *cn_Amat,
    double bc_0, double bc_1)



bool
update_solution_upwind15(int n, double *curr, double const *last,
    double alpha, double dx, double dt,
    double bc_0, double bc_1)
void 
compute_exact_solution(int n, double *a, double dx, char const *ic,
    double alpha, double t, double bc0, double bc1)
bool                          // false if unstable, true otherwise
update_solution_ftcs(
    int n,                    // number of samples
    double *uk1,              // new array of u(x,k+1) to compute/return
    double const *uk0,        // old/last array u(x,k) of samples computed
    double alpha,             // thermal diffusivity
    double dx, double dt,     // spacing in space, x, and time, t.
    double bc0, double bc1)   // boundary conditions @ x=0 & x=Lx


static void
initialize(void)

int finalize(int ti, double maxt, double change)

static bool
update_solution()

static double
update_output_files(int ti)

