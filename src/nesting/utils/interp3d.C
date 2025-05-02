#include "nesting_utils.h"

using namespace std;

void interp3d
(
  const std::vector<double>& x0, 
  const std::vector<double>& y0, 
  const std::vector<double>& z0,
  const std::vector<double>& u0,
  const double x1,
  const double y1,
  const double z1
)
{
  // Assuming this is a structured grid, monotonically incrasing in value
  // The problem with this grid is that it is not rectilinear,
  // therefore it is rather important to first construct a grid and then do the
  // interpolation.
  // The grid generation can be done 
}
