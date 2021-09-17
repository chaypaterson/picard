#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* This program tests an algorithm for a numerical approximation to Picard
 * iteration (with Euler as an initial guess) to solve a first-order ordinary
 * differential equation and initial value problem of the form
 *      dy/dx = f(x,y)
 * Instead of caching the entire function y(x) at many values of x to perform
 * the integration, we only store a local state that is just enough to move
 * foward by a given dx.
 * The algorithm looks like:
 *      1. Make a new Euler estimate using the last value, the function f, and
 *         the interval dx.
 *      2. Use the last estimate of y and current estimate to estimate
 *         \int f(x,y) dx over the interval. We can then update the first Picard
 *         estimate.
 *      3. Use this new estimate of y (and the last estimate of y in the hierarchy) 
 *         to update the new estimate of y.
 *      4. Repeat until convergence.
 *      5. Update x --> x + dx.
 *
 * We don't need to store the whole history of
 * Picard estimates. We can just use the best Picard estimate from last time for
 * the next Euler estimate. Then we only need to store:
 *      (x, last_best_estimate_of_y_at_x)
 *
 * This being said, it should become obvious that the algorithm is just a
 * generalised predictor-corrector method, related to Heun's method.
 * 
 * To test this algorithm, we will try to compute e. Adjusting dx provides a
 * nice demonstration of the tradeoff between truncation and roundoff errors.
 */

typedef double real; // FIXME this was intended to provide quasi-polymorphism
typedef real (*rhs_func)(real x, real y); // 
typedef void method; // methods should only have side effects and return nothing

typedef struct {
    // geometrically, the "state" of the ODE looks like a point in space: which
    // makes sense, because ODEs are basically vector fields. We are integrating
    // a modified ODE of the form
    //     dy/dt = f(x,y)
    //     dx/dt = 1
    // This naturally raises the question of how to deal with singular points
    // and singular forms of f like f = y^2. TODO.
    real x;
    real y;
} iter_state;

method one_step (real dx, rhs_func f, iter_state* state) {
    real ylast = state->y; // store the last value

    state->y += dx * f(state->x, ylast); // perform one Euler step

    // now iterate Picard until the new value converges:
    real pylast = ylast;
    while (pylast != state->y) {
        // one Picard step:
        pylast = state->y;
        // note: simpson's rule is even better here TODO
        state->y = ylast + 0.5 * (
                            f(state->x, ylast) + 
                            f(state->x + dx, state->y)
                        ) * dx;
    }

    state->x += dx;
}

real f(real x, real y) {
	// the rhs function we are integrating
	return y;
}

int main(void) {
    // set initial values
    iter_state xy;
    xy.x = 0.0;
    xy.y = 1.0;

    // choose a dx;
    real dx = 1e-6;

    // final value:
    real xfinal = 1.0;

    // I love one line while loops because I have been reading too much K&R:
    while (xy.x < xfinal - dx) one_step(dx, f, &xy);
    printf("%1.20f %1.20f\n", xy.x, xy.y);

    // correct last prediction for overshoot:
    real x_corr = xfinal - xy.x;
    one_step(x_corr, f, &xy);
    printf("%1.20f %1.20f\n", xy.x, xy.y);

    // compare value from math.h:
    printf("%1.20f %1.20f\n", xfinal, exp(xfinal));

    return 0;
}
