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
            // Define the callback function for this equation
            static void funcToSolve1(int n, IntPtr x, IntPtr fx)      // This is the function signature
            {
                double[] x1 = Fsolve.ExtractArrayFromPointer(n, x);   // Make an array for 'x' values from its Pointer
                double[] fx1 = Fsolve.ExtractArrayFromPointer(n, fx); // Make an array for 'fx' equation/function values from its Pointer
                fx1[0] = x1[0] - 1.0;                                 // Write equations/functions as f(x) = 0, here it becomes x - 1 = 0
                Fsolve.CopyArrayToPointer(n, fx1, fx);                // Copy fx1 array values to fx Pointer
            }

            // STEP 2:
            // Solve the function
            Fsolve.FunctionToSolve func1 = new (funcToSolve1);        // Wrap function so it can be called
            int unknownVariables1 = 1;                                // Give number of variables
            double[] xGuess1 = { -1.0 };                              // Give a guess value
            (double[] soln1, double[] fx1, string solutionCode1) = Fsolve.Fsolver(func1, unknownVariables1, xGuess1, Tolerance); // Call solver
            // Returns solution:
            // soln1 = Array containing solution
            // fx1 = Values of equations at soln1 (should be close to zero within Tolerance)
            // solutionCode1 = String providing information on exit code
            Fsolve.PrintArray ("xSolution", soln1, 8);                // Prints the solution to '8' decimals



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
            // Define the callback function for this system
            static void funcToSolve2(int n, IntPtr x, IntPtr fx)      // This is the function signature
            {
                double[] x1 = Fsolve.ExtractArrayFromPointer(n, x);   // Make an array for 'x' values from its Pointer
                double[] fx1 = Fsolve.ExtractArrayFromPointer(n, fx); // Make an array for 'fx' equation/function values from its Pointer
                fx1[0] = -x1[0] + 3.0 * x1[1] + 7.0 * x1[2];          // Write equations/functions as f(x) = 0
                fx1[1] = 2.0*x1[0] - 2.0*x1[1] - x1[2];
                fx1[2] = x1[0] + x1[1] + x1[2] - 1.0;
                Fsolve.CopyArrayToPointer(n, fx1, fx);                // Copy fx1 array values to fx Pointer
            }

            // STEP 2:
            // Solve the function
            Fsolve.FunctionToSolve func2 = new(funcToSolve2);        // Wrap function so it can be called
            int unknownVariables2 = 3;                               // Give number of variables
            double[] xGuess2 = { 0.0, 0.0, 0.0};                     // Give a guess value
            (double[] soln2, double[] fx2, string solutionCode2) = Fsolve.Fsolver(func2, unknownVariables2, xGuess2, Tolerance); // Call solver
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
            // Define the callback function for this system
            static void funcToSolve3(int n, IntPtr x, IntPtr fx)      // This is the function signature
            {
                double[] x1 = Fsolve.ExtractArrayFromPointer(n, x);   // Make an array for 'x' values from its Pointer
                double[] fx1 = Fsolve.ExtractArrayFromPointer(n, fx); // Make an array for 'fx' equation/function values from its Pointer
                fx1[0] = x1[0] + x1[1] - 1.0;                         // Write equations/functions as f(x) = 0
                fx1[1] = x1[0] * x1[0] - x1[1] - 5.0;
                Fsolve.CopyArrayToPointer(n, fx1, fx);                // Copy fx1 array values to fx Pointer
            }

            // STEP 2:
            // Solve the function
            Fsolve.FunctionToSolve func3 = new(funcToSolve3);        // Wrap function so it can be called
            int unknownVariables3 = 2;                               // Give number of variables
            double[] xGuess3a = { 0.0, 0.0 };                        // Give a guess value
            // Two separate guess values will be required to obtain both the solutions.
            // This guess will give the solution x = 2, y = -1
            (double[] soln3a, double[] fx3a, string solutionCode3a) = Fsolve.Fsolver(func3, unknownVariables3, xGuess3a, Tolerance); // Call solver
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
            (double[] soln3b, double[] fx3b, string solutionCode3b) = Fsolve.Fsolver(func3, unknownVariables3, xGuess3b, Tolerance); // Call solver
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
            // Define the callback function for this system
            static void funcToSolve4(int n, IntPtr x, IntPtr fx)      // This is the function signature
            {
                double[] x1 = Fsolve.ExtractArrayFromPointer(n, x);   // Make an array for 'x' values from its Pointer
                double[] fx1 = Fsolve.ExtractArrayFromPointer(n, fx); // Make an array for 'fx' equation/function values from its Pointer
                for (int k = 0; k < n; k++)                           // Write equations/functions as f(x) = 0
                {
                    double temp = (3.0 - 2.0 * x1[k]) * x1[k];
                    double temp1 = k != 0 ? x1[k - 1] : 0.0;
                    double temp2 = k != n - 1 ? x1[k + 1] : 0.0;
                    fx1[k] = temp - temp1 - 2.0 * temp2 + 1.0;
                }
                Fsolve.CopyArrayToPointer(n, fx1, fx);                // Copy fx1 array values to fx Pointer
            }
            // STEP 2:
            // Solve the function
            Fsolve.FunctionToSolve func4 = new(funcToSolve4);        // Wrap function so it can be called
            int unknownVariables4 = 9;                               // Give number of variables
            double[] xGuess4 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };  // Give a guess value
            (double[] soln4, double[] fx4, string solutionCode4) = Fsolve.Fsolver(func4, unknownVariables4, xGuess4, Tolerance); // Call solver
            // Returns solution:
            // soln4 = Array containing solution
            // fx4 = Values of equations at soln4 (should be close to zero within Tolerance)
            // solutionCode4 = String providing information on exit code
            Fsolve.PrintArray("xSolution", soln4, 7);                // Prints the solution to '7' decimals
        }
    }
}