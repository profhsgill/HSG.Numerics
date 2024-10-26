# HSG.Numerics

HSG.Numerics provides a solver for a system of nonlinear and linear equations with 'n' unknowns. It is inspired by the ***fsolve*** function of SciPy and MATLAB, and can be called from all .Net languages: F Sharp(F#), C Sharp (C#), and Visual Basic (VB).


## Motivation and Background

- SciPy/Python and MATLAB have a function called ***fsolve*** to solve a sytem of nonlinear equations. It is based on [MINPACK](https://www.mcs.anl.gov/~more/ANL8074a.pdf).

- HSG.Numerics provides a similar ***fsolve*** function using the C implementation of [fsolve](https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html) by John Burkardt, which is also based on MINPACK. 

- This function can be called by all the .Net languages: F#, C#, and VB.


## Version notice

Starting version 1.0.5 the usage has been simplified. The user has to just define the function to be solved and then pass it to the solver. Other aspects of handling C interoperability are handled by the library.


## Installation

Use nuget to install.

```bash
dotnet add package HSG.Numerics
```

`OR` 

```bash
In Visual Studio right click on project in Solution Explorer and in the menu that appears click on 'Manage NuGet Packages'. Then 'browse' to find HSG.Numerics and install it.
```



## Compatibility 

- The dlls were compiled for .Net 6 on Windows 10. But they should likely work on Windows 11. 

- Library has not been checked for other Operating Systems.

- fsharp.core with version number at least 8.0.100 is required (because dlls were compiled with this version).



## Known issues

When the library is included in an interactive script file for F# as follows:

 `#r "nuget: HSG.Numerics, 1.0.5"` (any version)
 
 an error gets raised that reads "bad cli header, rva 0".

 The library must be added into a .fs file in a "non interactive" manner and compiled.




## Examples

One example is given below for F#, C# and VB. More can be found at project github repo [HSG.Numerics](https://github.com/profhsgill/HSG.Numerics)


`F#`

```F#
open System
open HSG.Numerics
[<EntryPoint>]
let main argv =
    let Tolerance = 1e-10 // tolerance for all test cases
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
    let funcToSolve4 (x: double[]) =  // This is the function signature          
        let n = x.Length
        let fx = Array.zeroCreate n
        for k = 0 to n-1 do                              // Write equations/functions as f(x) = 0
            let temp = (3.0 - 2.0*x[k])*x[k]
            let temp1 = if k <> 0 then x[k - 1] else 0.0
            let temp2 = if k <> n-1 then x[k + 1] else 0.0
            fx[k] <- temp - temp1 - 2.0*temp2 + 1.0
        fx // Return fx array
    
    // STEP 2:
    // Solve the function
    let unknownVariables4 = 9                              // Give number of variables 
    let xGuess4:double array = Array.zeroCreate 9          // Give a guess value
    let solveResult4 = Fsolve.Fsolver(funcToSolve4, unknownVariables4, xGuess4, Tolerance) // Call solver
    let (soln4, fx4, solutionCode4) = solveResult4        // Returns solution:
    // soln4 = Array containing solution
    // fx4 = Values of equations at soln4 (should be close to zero within Tolerance)
    // solutionCode4 = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln4 7                  // Prints the solution to '7' decimals

    0

```


`C#`

```C#
using HSG.Numerics;

namespace FsolveTester
{
    class Tester
    {
        static void Main(string[] args)
        {
            double Tolerance = 1e-10; // tolerance for all test cases

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

        }
    }
}

```


`VB`

```VB
Imports HSG.Numerics
Module FsolveTester

    Sub Main()

        Dim Tolerance As Double = 0.0000000001 ' Tolerance for all test cases

       
        ' Solve function 4
        ' ----------------
        Dim unknownVariables4 As Integer = 9                            ' Give number of variables 
        Dim xGuess4() = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}   ' Give a guess value
        Dim solveResult4 = Fsolve.Fsolver(AddressOf funcToSolve4, unknownVariables4, xGuess4, Tolerance) ' Call solver
        ' Returns solution:
        ' soln4 = Array containing solution
        ' fx4 = Values of equations at soln4 (should be close to zero within Tolerance)
        ' solutionCode4 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult4.Item1, 7)

    End Sub
    
    
    ' ============================================
    ' ....... TEST 4 .............................
    ' ....... Non Linear System ..................
    ' ....... Unknown variables = 9 ..............
    ' ============================================
    '
    ' This is an example from original MINPACK User Guide
    '
    ' Equations:
    '
    ' This is a tri-diagonal matrix
    '
    ' (3 – 2*x(0)) * x(0)                     -2*x(1)                             = -1
    '             -x(i-1)      +      (3-2*x(i))*x(i)                   -2*x(i+1) = -1, i=1,7
    '                                           -x(7)      +      (3-2*x(8))*x(8) = -1
    '
    ' Original solution from User Guide:
    ' 
    ' Solutions:
    ' -0.5706545      -0.6816283      -0.7017325
    ' -0.7042129      -0.7013690      -0.6918656
    ' -0.6657920      -0.5960342      -0.4164121
    Function funcToSolve4(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        For k As Integer = 0 To n - 1                                   ' Write equations/functions as f(x) = 0
            Dim temp As Double = (3.0 - 2.0 * x(k)) * x(k)
            Dim temp1 As Double = If((k <> 0), x(k - 1), 0.0)
            Dim temp2 As Double = If((k <> n - 1), x(k + 1), 0.0)
            fx(k) = temp - temp1 - 2.0 * temp2 + 1.0
        Next
        Return fx
    End Function

End Module

```


## License

[MIT](https://choosealicense.com/licenses/mit/)

## References

1. [MINPACK-1](https://www.mcs.anl.gov/~more/ANL8074a.pdf)
2. [John Burkardt - fsolve](https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html)
3. .NET 2.0 Interoperability Recipes: A Problem-Solution Approach (Expert's Voice in .NET) by  Bruce Bukovics ( ISBN-13 : 978-1590596692 )
