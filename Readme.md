# HSG.Numerics

HSG.Numerics provides a solver for a system of nonlinear and linear equations with 'n' unknowns. It is inspired by the ***fsolve*** function of SciPy and MATLAB, and can be called from all .Net languages (F#, C#, VB).


## Motivation and Background

- SciPy/Python and MATLAB have a function called ***fsolve*** to solve a sytem of nonlinear equations. It is based on [MINPACK](https://www.mcs.anl.gov/~more/ANL8074a.pdf).

- HSG.Numerics provides a similar ***fsolve*** function using the C implementation of [fsolve](https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html), which is also based on MINPACK. 

- This function can be called by all the .Net languages: F#, C#, and VB.


## Version notice

Starting version 1.0.4 and above the names of the following functions have been changed.

- *Fsolve.ExtractArrayFromPointer* is now called *Fsolve.MakeArray*
- *Fsolve.CopyArrayToPointer* is now called *Fsolve.CopyArray*


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

  `#r "nuget: HSG.Numerics, 1.0.3"` (any version)
 
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
    // Define the callback function for this system
    let funcToSolve4 (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature
        let x1 = Fsolve.MakeArray n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx    // Make an array for 'fx' equation/function values from its Pointer                   
        for k = 0 to n-1 do                              // Write equations/functions as f(x) = 0
            let temp = (3.0 - 2.0*x1[k])*x1[k]
            let temp1 = if k <> 0 then x1[k - 1] else 0.0
            let temp2 = if k <> n-1 then x1[k + 1] else 0.0
            fx1[k] <- temp - temp1 - 2.0*temp2 + 1.0
        Fsolve.CopyArray n fx1 fx               // Copy fx array values to fx Pointer
        ()                                      // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func4 = Fsolve.FunctionToSolve (funcToSolve4)      // Wrap function so it can be called
    let unknownVariables4 = 9                              // Give number of variables 
    let xGuess4:double array = Array.zeroCreate 9          // Give a guess value
    let solveResult4 = Fsolve.Fsolver(func4, unknownVariables4, xGuess4, Tolerance) // Call solver
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
            // Define the callback function for this system
            static void funcToSolve4(int n, IntPtr x, IntPtr fx)      // This is the function signature
            {
                double[] x1 = Fsolve.MakeArray(n, x);   // Make an array for 'x' values from its Pointer
                double[] fx1 = Fsolve.MakeArray(n, fx); // Make an array for 'fx' equation/function values from its Pointer
                for (int k = 0; k < n; k++)                           // Write equations/functions as f(x) = 0
                {
                    double temp = (3.0 - 2.0 * x1[k]) * x1[k];
                    double temp1 = k != 0 ? x1[k - 1] : 0.0;
                    double temp2 = k != n - 1 ? x1[k + 1] : 0.0;
                    fx1[k] = temp - temp1 - 2.0 * temp2 + 1.0;
                }
                Fsolve.CopyArray(n, fx1, fx);                // Copy fx1 array values to fx Pointer
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

```


`VB`

```VB
Imports HSG.Numerics
Module FsolveTester

    Sub Main()

        Dim Tolerance As Double = 0.0000000001 ' Tolerance for all test cases

       
        ' Solve function 4
        ' ----------------
        Dim func4 As New Fsolve.FunctionToSolve(AddressOf funcToSolve4) ' Wrap function so it can be called
        Dim unknownVariables4 As Integer = 9                            ' Give number of variables 
        Dim xGuess4() = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}   ' Give a guess value
        Dim solveResult4 = Fsolve.Fsolver(func4, unknownVariables4, xGuess4, Tolerance) ' Call solver
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
    Sub funcToSolve4(n As Integer, x As IntPtr, fx As IntPtr)           ' This is the function signature
        Dim x1() As Double = Fsolve.MakeArray(n, x)       ' Make an array for 'x' values from its Pointer
        Dim fx1() As Double = Fsolve.MakeArray(n, fx)     ' Make an array for 'fx' equation/function values from its Pointer
        For k As Integer = 0 To n - 1                                   ' Write equations/functions as f(x) = 0
            Dim temp As Double = (3.0 - 2.0 * x1(k)) * x1(k)
            Dim temp1 As Double = If((k <> 0), x1(k - 1), 0.0)
            Dim temp2 As Double = If((k <> n - 1), x1(k + 1), 0.0)
            fx1(k) = temp - temp1 - 2.0 * temp2 + 1.0
        Next
        Fsolve.CopyArray(n, fx1, fx)                           ' Copy fx array values to fx Pointer
    End Sub

End Module

```


## License

[MIT](https://choosealicense.com/licenses/mit/)

## References

1. [MINPACK-1](https://www.mcs.anl.gov/~more/ANL8074a.pdf)
2. [John Burkardt - fsolve](https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html)
3. .NET 2.0 Interoperability Recipes: A Problem-Solution Approach (Expert's Voice in .NET) by  Bruce Bukovics ( ISBN-13 : 978-1590596692 )
