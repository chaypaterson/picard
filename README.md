# picard

This project is an implementation of a numerical integrator in C.

The numerical integration routine is a demonstration of both:
    1. Numerical integration
    2. Picard iteration

The integrator acts on an initial value problem of the form
    dy/dx = f(x,y(x))
    y(x = initial_x_value) = initial_y_value
and produces the best possible approximation of y(x == chosen_final_value) given
a finite step size dx.

This "best possible approximation" is arrived at by iterating
    y(x) = y_0 + \int_{x' = 0}^x f(x',y(x')) dx'
or equivalently:
    y(x + \Delta x) = y(x) + \int_{z = 0}^{\Delta x} f(x + z, y(x + z)) dz
until it converges (to within the available precision). The latter integral is
what is approximated numerically: the precise method can be adjusted.

The fixed point of the numerical integration should be a good approximation to
the solution of the underlying differential equation. It is important that the
numerical integration rule should depend on both the last estimate "y(x)" and
the current estimate "y(x + \Delta x)" in order to converge.

I have tried versions of the midpoint rule and Simpson's rule with various step
sizes dx.
