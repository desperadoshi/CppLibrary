# Cubic Spline C++ Implementation

## Math

Refer to Numerical Recipe. The selected pages are attached. [PDF link](./Cubic_Spline_Interpolation.pdf)

The "Not-A-Knot" boundary condition is not explained in the above PDF file. Refer this [note](./not-a-knot_bc.pdf). In fact, the "Not-A-Knot" BC is to set the 3rd derivative to be continuous at the 2nd point and also the 2nd last point.

In fact, there is another point of view of thinking about the formulation of the cubic spline. The way in the Numerical Recipe is perfect, while another way is from how many equations to match the number of unknowns. Such method leads difficulties of applying the "Not-A-Knot" BC. Because the unknowns are not the 2nd derivatives, which is not straightforward from the viewpoint of the derivatives.

## C++

An OOP implementation.

The procedure is simple.

-   Assemble the coefficients matrix $A$ and the RHS $b$ to form $A x = b$
-   Solve the system equations to get $x$, i.e. the 2nd derivative at each point
-   Given a target $x$, apply the overloaded operator `()` to obtain the target interpolation value.

## CMake

`CMakeLists.txt` is organized as follows,

-   Set the project name
-   Set the CMake module path
    -   We use a cmake file to find the external library. We need to tell the master `CMakeLists.txt` where the cmake module is located.
-   Find the external library by the cmake module file
    -   It is better to use the cmake module shipped within the external library.
-   Include all the codes
-   Set necessary compiler flags
-   Assemble the final binary file

It is also possible to organize in another way. Compile the `CubicSpline` class to have a static library, which could be used as an external library. Currently, I choose the 1st way.

## Eigen3

The Eigen3 package is employed to solve the linear system. Only the function `SolveCoef` contains codes from the Eigen3. If it needs external library to solve the linear system, it is easy to change.

The Eigen3 package is easy to use. Include the header file, most commonly `<Eigen/Dense>`. The matrix is allocated by `Eigen::MatrixXd A(n,n);` and the vector is given by `Eigen::VectorXd b(n);`. To solve the linear system of equations, the class member function is applied, `X=A.colPivHouseholderQr().solve(b);`.

## Source Codes

The source codes are hosted in [Github](https://github.com/desperadoshi/CubicSpline).

