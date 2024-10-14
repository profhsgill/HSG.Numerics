using HSG.Numerics;


namespace FsolveTester
{
    class Tester
    {
        static void Main(string[] args)
        {
            double Tolerance = 1e-10; // tolerance for all test cases

            // DEFINE TESTS


            // ============================================
            // ....... TEST 1 .............................
            // ....... Linear System ......................
            // ....... Unknown variables = 1 ..............
            // ============================================
            //
            // Equations:
            // 1) x = 1
            //
            // Solution:
            // x = 1

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve1(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                fx[0] = x[0] - 1.0;   // Write equations/functions as f(x) = 0, here it becomes x - 1 = 0
                return fx;
            }

            // STEP 2:
            // Solve the function
            int unknownVariables1 = 1;                                // Give number of variables
            double[] xGuess1 = { -1.0 };                              // Give a guess value
            (double[] soln1, double[] fx, string solutionCode1) = Fsolve.Fsolver(funcToSolve1, unknownVariables1, xGuess1, Tolerance); // Call solver
            // Returns solution:
            // soln1 = Array containing solution
            // fx = Values of equations at soln1 (should be close to zero within Tolerance)
            // solutionCode1 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln1, 8);                // Prints the solution to '8' decimals



            // ============================================
            // ....... TEST 2 .............................
            // ....... Linear System ......................
            // ....... Unknown variables = 3 ..............
            // ============================================
            //
            // Equations:
            // -x + 3y + 7z = 0
            // 2x - 2y - z = 0
            // x + y + z = 1
            //
            // Solution:
            // x = 0.55, y = 0.65, and z = -0.2

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve2(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                fx[0] = -x[0] + 3.0 * x[1] + 7.0 * x[2];          // Write equations/functions as f(x) = 0
                fx[1] = 2.0 * x[0] - 2.0 * x[1] - x[2];
                fx[2] = x[0] + x[1] + x[2] - 1.0;
                return fx;
            }

            // STEP 2:
            // Solve the function
            int unknownVariables2 = 3;                               // Give number of variables
            double[] xGuess2 = { 0.0, 0.0, 0.0 };                     // Give a guess value
            (double[] soln2, double[] fx2, string solutionCode2) = Fsolve.Fsolver(funcToSolve2, unknownVariables2, xGuess2, Tolerance); // Call solver
            // Returns solution:
            // soln2 = Array containing solution
            // fx2 = Values of equations at soln2 (should be close to zero within Tolerance)
            // solutionCode2 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln2, 8);                // Prints the solution to '8' decimals



            // ============================================
            // ....... TEST 3 .............................
            // ....... Non Linear System ..................
            // ....... Unknown variables = 2 ..............
            // ============================================
            //
            // Equations:
            // x + y = 1
            // y = x**2 - 5
            //
            // Solutions:
            // (1) x = -3, y = 4;      (2) x = 2, y = -1

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve3(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                fx[0] = x[0] + x[1] - 1.0;                         // Write equations/functions as f(x) = 0
                fx[1] = x[0] * x[0] - x[1] - 5.0;
                return fx;
            }

            // STEP 2:
            // Solve the function
            int unknownVariables3 = 2;                               // Give number of variables
            double[] xGuess3a = { 0.0, 0.0 };                        // Give a guess value
            // Two separate guess values will be required to obtain both the solutions.
            // This guess will give the solution x = 2, y = -1
            (double[] soln3a, double[] fx3a, string solutionCode3a) = Fsolve.Fsolver(funcToSolve3, unknownVariables3, xGuess3a, Tolerance); // Call solver
            // Returns solution:
            // soln3a = Array containing solution
            // fx3a = Values of equations at soln3a (should be close to zero within Tolerance)
            // solutionCode3a = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln3a, 8);               // Prints the solution to '8' decimals
            // Give 2nd guess
            //----------------
            double[] xGuess3b = { -2.0, 0.0 };                       // Give a guess value
            // Two separate guess values will be required to obtain both the solutions.
            // This guess will give the solution x = -3, y = 4
            (double[] soln3b, double[] fx3b, string solutionCode3b) = Fsolve.Fsolver(funcToSolve3, unknownVariables3, xGuess3b, Tolerance); // Call solver
            // Returns solution:
            // soln3b = Array containing solution
            // fx3b = Values of equations at soln3b (should be close to zero within Tolerance)
            // solutionCode3b = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln3b, 8);               // Prints the solution to '8' decimals



            // ============================================
            // ....... TEST 4 .............................
            // ....... Non Linear System ..................
            // ....... Unknown variables = 9 ..............
            // ============================================
            //
            // This is an example from original MINPACK User Guide
            //
            // Equations:
            //
            // This is a tri-diagonal matrix
            //
            // (3 – 2*x(0)) * x(0)                     -2*x(1)                             = -1
            //             -x(i-1)      +      (3-2*x(i))*x(i)                   -2*x(i+1) = -1, i=1,7
            //                                           -x(7)      +      (3-2*x(8))*x(8) = -1
            //
            // Original solution from User Guide:
            // 
            // Solutions:
            // -0.5706545      -0.6816283      -0.7017325
            // -0.7042129      -0.7013690      -0.6918656
            // -0.6657920      -0.5960342      -0.4164121

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve4(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                for (int k = 0; k < n; k++)                           // Write equations/functions as f(x) = 0
                {
                    double temp = (3.0 - 2.0 * x[k]) * x[k];
                    double temp1 = k != 0 ? x[k - 1] : 0.0;
                    double temp2 = k != n - 1 ? x[k + 1] : 0.0;
                    fx[k] = temp - temp1 - 2.0 * temp2 + 1.0;
                }
                return fx;
            }
            // STEP 2:
            // Solve the function
            int unknownVariables4 = 9;                               // Give number of variables
            double[] xGuess4 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };  // Give a guess value
            (double[] soln4, double[] fx4, string solutionCode4) = Fsolve.Fsolver(funcToSolve4, unknownVariables4, xGuess4, Tolerance); // Call solver
            // Returns solution:
            // soln4 = Array containing solution
            // fx4 = Values of equations at soln4 (should be close to zero within Tolerance)
            // solutionCode4 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln4, 7);                // Prints the solution to '7' decimals



            // ============================================
            // ....... TEST 5 .............................
            // ....... Non Linear System ..................
            // ....... Unknown variables = 2 ..............
            // ============================================
            //
            // This is an example from https://www.mathworks.com/help/optim/ug/fsolve.html
            //
            // Equations:
            //
            // e**(−e**(−(x+x2))) = x2*(1+x**2)
            // x*cos(x2) + x2*sin(x) = 0.5
            // 
            // Solutions:
            // 0.3532    0.6061

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve5(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                fx[0] = Math.Exp(-Math.Exp(-(x[0] + x[1]))) - x[1] * (1.0 + x[0] * x[0]); // Write equations/functions as f(x) = 0
                fx[1] = x[0] * Math.Cos(x[1]) + x[1] * Math.Sin(x[0]) - 0.5;
                return fx;
            }
            // STEP 2:
            // Solve the function
            int unknownVariables5 = 2;                               // Give number of variables
            double[] xGuess5 = { 0.0, 0.0 };  // Give a guess value
            (double[] soln5, double[] fx5, string solutionCode5) = Fsolve.Fsolver(funcToSolve5, unknownVariables5, xGuess5, Tolerance); // Call solver
            // Returns solution:
            // soln5= Array containing solution
            // fx5 = Values of equations at soln5 (should be close to zero within Tolerance)
            // solutionCode5 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln5, 4);                // Prints the solution to '7' decimals



            // ============================================
            // ....... TEST 6 .............................
            // ....... Non Linear System ..................
            // ....... Unknown variables = 2 ..............
            // ============================================
            //
            // This is an example from https://www.mathworks.com/help/optim/ug/fsolve.html
            //
            // Equations:
            //
            // 2x1 − x2  = e**(−x)
            // −x + 2x2 = e**(−x2)
            //
            // Solutions:
            // 0.5671    0.5671

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve6(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                fx[0] = 2.0 * x[0] - x[1] - Math.Exp(-x[0]);   // Write equations/functions as f(x) = 0
                fx[1] = -x[0] + 2.0 * x[1] - Math.Exp(-x[1]);
                return fx;
            }
            // STEP 2:
            // Solve the function
            int unknownVariables6 = 2;                               // Give number of variables
            double[] xGuess6 = { 0.0, 0.0 };  // Give a guess value
            (double[] soln6, double[] fx6, string solutionCode6) = Fsolve.Fsolver(funcToSolve6, unknownVariables6, xGuess6, Tolerance); // Call solver
            // Returns solution:
            // soln6= Array containing solution
            // fx6 = Values of equations at soln6 (should be close to zero within Tolerance)
            // solutionCode6 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln6, 4);                // Prints the solution to '7' decimals



            // ============================================
            // ....... TEST 7 .............................
            // ....... Non Linear System ..................
            // ....... Unknown variables = 1 ..............
            // ============================================
            //
            // This is an example from https://www.mathworks.com/help/optim/ug/fsolve.html
            //
            // Equations:
            //
            //  a(1+cos(b))**2 = 4x * exp(-2c*(a-x)**2) 
            //  
            //
            // Solutions:
            // Depends on a, b, c
            // 0.70548923 for a = b = c = 1

            // STEP 1: 
            // Define the function for this equation
            static double[] funcToSolve7(double[] x)      // This is the function signature
            {
                int n = x.Length;
                double[] fx = new double[n];
                double a = 1.0;
                double b = 1.0;
                double c = 1.0;
                // Write equations/functions as f(x) = 0
                double term1 = 1.0 + Math.Cos(b);
                double term2 = -2.0 * c * (a - x[0]) * (a - x[0]);
                fx[0] = a * term1 * term1 - 4.0 * x[0] * Math.Exp(term2);
                return fx;
            }
            // STEP 2:
            // Solve the function
            int unknownVariables7 = 1;                               // Give number of variables
            double[] xGuess7 = { 0.0 };  // Give a guess value
            (double[] soln7, double[] fx7, string solutionCode7) = Fsolve.Fsolver(funcToSolve7, unknownVariables7, xGuess7, Tolerance); // Call solver
            // Returns solution:
            // soln7= Array containing solution
            // fx7 = Values of equations at soln6 (should be close to zero within Tolerance)
            // solutionCode7 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln7, 8);                // Prints the solution to '7' decimals
        }
    }
}