﻿// <copyright file="HSG.Numerics.fs">
//
// Copyright (c) 2024 Harvinder Singh Gill <profhsgill@gmail.com>
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>


// -------
// Authors
// -------

// Original FORTRAN77 version of MINPACK written by:
// Jorge More, Danny Sorenson, Burton Garbow, Kenneth Hillstrom.

// The C version was written by John Burkardt under GNU LGPL License.
// https://people.sc.fsu.edu/~jburkardt/c_src/fsolve/fsolve.html

// The following code is written by Harvinder Singh Gill.
// The following code dynamically links to C code by John Burkardt.

// -------
// Date
// -------
// November 12, 2023

// ------------
// Revised Date
// ------------
// October 13, 2024

// REVISION
// --------
// Improved the library such that the user ONLY provides function to be solved
// The function is written in fsharp, csharp or VB



namespace HSG.Numerics

open System
open System.Runtime.InteropServices


/// <summary>
/// Library for .Net to solve a system of non-linear equations (N equations in N unknowns).
/// Inspired by fsolve() function in SciPY and MATLAB, and based on MINPACK package written originally in FORTRAN.
/// </summary>
module Fsolve =

    /// <summary>
    /// Delegate of function to be solved given by user.
    /// </summary>
    type UserFunction = delegate of double[] -> double[]

    /// <summary>
    /// Delegate of function required by C fsolve function.
    /// </summary>
    type CNativeFunction = delegate of int * IntPtr * IntPtr -> unit
    


    [<DllImport(@"C.Numerics.dll",
    CallingConvention=CallingConvention.Cdecl)>]
    extern int fsolve(CNativeFunction CfuncToSolve, int unknownVariableCount, double[] x, double[] fx, double tolerance)

    // HELPER FUNCTIONS
    // ----------------
    // ----------------
    /// <summary>
    /// Given a Pointer (address), get the array values
    /// </summary>
    /// <param name="arrayElementCount"> Number of array elements </param>
    /// <param name="pointer"> Pointer to starting location of the Array</param>
    /// <returns> Array of Double [] </returns>
    let MakeArray (arrayElementCount:int) (pointer:IntPtr) =
        let mutable newArray:double[] = Array.zeroCreate arrayElementCount
        Marshal.Copy(pointer, newArray, 0, arrayElementCount) // Copy pointer values to array
        newArray // Return array

    /// <summary>
    /// Given an Array, copy its values to a given Pointer (address)
    /// </summary>
    /// <param name="arrayElementCount"> Number of array elements </param>
    /// <param name="sourceArray"> Array containing values to be copied</param>
    /// <param name="destinationPointer"> Pointer(address) where Array is to be copied</param>
    /// <returns> No return </returns>
    let CopyArray (arrayElementCount:int) (sourceArray:double[]) (destinationPointer:IntPtr) =
        Marshal.Copy(sourceArray, 0, destinationPointer, arrayElementCount)

    /// <summary>
    /// Prints array values
    /// </summary>
    /// <param name="name"> Name of the array</param>
    /// <param name="arrayToPrint"> Array whose values are to be printed</param>
    /// <param name="decimalPlaces"> Decimal places to include in the print</param>
    /// <returns> No return </returns>
    let PrintArray (name:string) (arrayToPrint:double[]) (decimalPlaces:int) =
        let n = arrayToPrint.Length
        for i in [0..n-1] do
            printfn "Value %s[%d] = %.*f" name i decimalPlaces arrayToPrint[i] // * gets replaced by decimalPlaces
        printfn "========= Print complete ========="


    /// <summary>
    /// Function to be called from dotnet language (C#, F#, VB.Net) to solve non-linear equations.
    /// </summary>
    /// <param name="func"> Function to be solved. It accepts an array of 'x' values and returns computed equation values at these values</param>
    /// <param name="unknownVariableCount"> Number of unknowns in the system </param>
    /// <param name="xGuess"> Array of initial guess values </param>
    /// <param name="tolerance"> Acceptable error between consecutive iterations of function values </param>
    /// <returns> Tuple of (xSolution, functionValuesAtSolution, solverReturnCode)</returns>
    // solverReturnCodes are translated into text as follows:
    // 0: Improper input parameters.
    // 1: Success: Relative error between consecutive iterations is less than or equal to tolerance.
    // 2: Number of calls to function has reached or exceeded 200*(n+1).
    // 3: Tolerance is too small. No further improvement in the approximate solution x is possible.
    // 4: Iteration is not making good progress.
    let Fsolver(func:UserFunction, unknownVariableCount:int, xGuess:double[], tolerance:double) =
        // func:UserFunction cannot be directly consumed for calling C fsolve function
        // The following converts 'func' into a form so it can be used for calling C fsolve function
        let CfuncToSolve (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature needed by 'C' fsolve functiion
            let x1 = MakeArray n x      // Make an array for 'x' values from its Pointer                 
            let fx1 = func.Invoke x1    // Call user function using this 'x' values array
            CopyArray n fx1 fx          // Copy the fx1 array values to fx Pointer
            ()      // returns unit mimicking the C implementation return of void from called function

        let mutable xSolution = Array.zeroCreate unknownVariableCount
        xSolution <- xGuess
        let mutable functionValuesAtSolution = Array.zeroCreate unknownVariableCount

        // C fsolve function changes xSolution and functionValuesAtSolution in every iteration
        // So they are made muatble

        let resultCode = fsolve(CfuncToSolve, unknownVariableCount, xSolution, functionValuesAtSolution, tolerance)
        let AnalyzeResult resultCode =
            if resultCode = 0 then "Improper input parameters."
            elif resultCode = 1 then "Success: Relative error between consecutive iterations is less than or equal to tolerance."
            elif resultCode = 2 then "Number of calls to function has reached or exceeded 200*(n+1), where n = number of unknowns."
            elif resultCode = 3 then "Tolerance is too small. No further improvement in the approximate solution x is possible."
            elif resultCode = 4 then "Iteration is not making good progress."
            else "Some unexpected error has occurred. Please check your function."
        xSolution, functionValuesAtSolution, AnalyzeResult resultCode
