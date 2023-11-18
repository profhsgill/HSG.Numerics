# HSG.Numerics

HSG.Numerics provides a solver for a system of nonlinear and linear equations with 'n' unknowns. It is inspired by the ***fsolve*** function of SciPy and MATLAB, and can be called from all .Net languages (F#, C#, VB).


## Motivation and Background

- While learning F# I had to solve a system of nonlinear equations. I searched a lot but could not find an implementation of a nonlinear solver for a system of equations for F# nor for C#. 

- SciPy/Python and MATLAB have a function called ***fsolve*** to solve a sytem of nonlinear equations. I then started researching into ways to write an fsolve function for F#.

- Since the ***fsolve*** function of SciPy and MATLAB uses MINPACK, I started looking into numerous implementations of [MINPACK](https://www.mcs.anl.gov/~more/ANL8074a.pdf). Because MINPACK is written in FORTRAN, its direct consumption in F# is difficult so I looked at C language implementations of MINPACK. However, the C implementations of MINPACK became difficult to consume in F# due to lack of examples or custom data structures implemented by these libraries.

- I then stumbled upon the C implementation of [fsolve](https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html) by John Burkardt. This implementation uses built in data types of C without introducing custom data types, and it seemed plausible that I might be able to consume it in F#.

- I had to start learning Interoperability between native C code and F# code. It was a steep learning curve, but in the end I succeeded in writing a library that uses the John Burkardt fsolve code written in C and allows its use in F#. I wrote the F# code in a way that it can also be consumed in C# and VB.

- The most difficult part was to determine how a function written in F# that takes in guess values of unknown variable (say 'x') and returns an array of equation values at this 'x' can be passed to C code where the fsolve function can repeatedly call it in an iterating loop.

- The HSG.Numerics is a result of this effort. It can now be called by all the .Net languages: F#, C#, and VB.

- I hope this will be useful to others.


## Installation

Use nuget to install.

```bash
dotnet add package HSG.Numerics
```
`OR` for F# script files add the following at the top of the script file.
```bash
#r "nuget: HSG.Numerics, 1.0.0"
```

## Compatibility

The dlls were compiled for .Net 6 on Windows 10. But they should likely work on Windows 11. I have not checked other Operating Systems.



## Examples

One example for each F#, C# and VB is given below.

More examples can be found in the 'examples' folder

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
        let x1 = Fsolve.ExtractArrayFromPointer n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.ExtractArrayFromPointer n fx    // Make an array for 'fx' equation/function values from its Pointer                   
        for k = 0 to n-1 do                              // Write equations/functions as f(x) = 0
            let temp = (3.0 - 2.0*x1[k])*x1[k]
            let temp1 = if k <> 0 then x1[k - 1] else 0.0
            let temp2 = if k <> n-1 then x1[k + 1] else 0.0
            fx1[k] <- temp - temp1 - 2.0*temp2 + 1.0
        Fsolve.CopyArrayToPointer n fx1 fx               // Copy fx array values to fx Pointer
        ()                                               // Returns equivalent of void in 'C'
    
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
        Dim x1() As Double = Fsolve.ExtractArrayFromPointer(n, x)       ' Make an array for 'x' values from its Pointer
        Dim fx1() As Double = Fsolve.ExtractArrayFromPointer(n, fx)     ' Make an array for 'fx' equation/function values from its Pointer
        For k As Integer = 0 To n - 1                                   ' Write equations/functions as f(x) = 0XXXXXXXXXXXXXXXXXXXXXXXXX ZXZX
            Dim temp As Double = (3.0 - 2.0 * x1(k)) * x1(k)
            Dim temp1 As Double = If((k <> 0), x1(k - 1), 0.0)
            Dim temp2 As Double = If((k <> n - 1), x1(k + 1), 0.0)
            fx1(k) = temp - temp1 - 2.0 * temp2 + 1.0
        Next
        Fsolve.CopyArrayToPointer(n, fx1, fx)                           ' Copy fx array values to fx PointerC
    End Sub

End Module

```


## License

[MIT](https://choosealicense.com/licenses/mit/)

## References

1. [MINPACK-1](https://www.mcs.anl.gov/~more/ANL8074a.pdf)
2. [John Burkardt - fsolve](https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html)
3. .NET 2.0 Interoperability Recipes: A Problem-Solution Approach (Expert's Voice in .NET) by  Bruce Bukovics ( ISBN-13 : 978-1590596692 )
