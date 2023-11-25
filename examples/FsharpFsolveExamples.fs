open System
open HSG.Numerics
[<EntryPoint>]
let main argv =
    let Tolerance = 1e-10 // tolerance for all test cases
    
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
    let funcToSolve1 (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature
        let x1 = Fsolve.MakeArray n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx    // Make an array for 'fx' equation/function values from its Pointer
        fx1[0] <- x1[0] - 1.0                            // Write equations/functions as f(x) = 0, here it becomes x - 1 = 0
        Fsolve.CopyArray n fx1 fx               // Copy fx1 array values to fx Pointer
        ()                                               // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func1 = Fsolve.FunctionToSolve (funcToSolve1)     // Wrap function so it can be called
    let unknownVariables1 = 1                             // Give number of variables 
    let xGuess1:double array = [| -2.5 |]                 // Give a guess value
    let solveResult1 = Fsolve.Fsolver(func1, unknownVariables1, xGuess1, Tolerance) // Call solver
    let (soln1, fx1, solutionCode1) = solveResult1        // Returns solution:
    // soln1 = Array containing solution
    // fx1 = Values of equations at soln1 (should be close to zero within Tolerance)
    // solutionCode1 = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln1 8                 // Prints the solution to '8' decimals



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
    let funcToSolve2 (n:int) (x: IntPtr) (fx: IntPtr) =   // This is the function signature
        let x1 = Fsolve.MakeArray n x       // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx     // Make an array for 'fx' equation/function values from its Pointer
        fx1[0] <- -1.0 * x1[0] + 3.0 * x1[1] + 7.0 * x1[2]// Write equations/functions as f(x) = 0
        fx1[1] <- 2.0 * x1[0] - 2.0 * x1[1] - 1.0 * x1[2]
        fx1[2] <- x1[0] + x1[1] + x1[2] - 1.0
        Fsolve.CopyArray n fx1 fx                // Copy fx1 array values to fx Pointer
        ()                                                // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func2 = Fsolve.FunctionToSolve (funcToSolve2)     // Wrap function so it can be called
    let unknownVariables2 = 3                             // Give number of variables 
    let xGuess2:double array = [| 0.0; 0.0; 0.0 |]        // Give a guess value
    let solveResult2 = Fsolve.Fsolver(func2, unknownVariables2, xGuess2, Tolerance) // Call solver
    let (soln2, fx2, solutionCode2) = solveResult2        // Returns solution:
    // soln2 = Array containing solution
    // fx2 = Values of equations at soln2 (should be close to zero within Tolerance)
    // solutionCode2 = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln2 8                 // Prints the solution to '8' decimals



    // ============================================
    // ....... TEST 3 .............................
    // ....... Non Linear System ......................
    // ....... Unknown variables = 2 ..............
    // ============================================
    //
    // Equations:
    // x + y = 1
    // y = x**2 - 5
    // solution: [x=-3, y=4] and [x=2, y=-1] 
    //
    // Solutions:
    // (1) x = -3, y = 4;      (2) x = 2, y = -1

    // STEP 1: 
    // Define the callback function for this system
    let funcToSolve3 (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature
        let x1 = Fsolve.MakeArray n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx    // Make an array for 'fx' equation/function values from its Pointer
        fx1[0] <- x1[0] + x1[1] - 1.0                    // Write equations/functions as f(x) = 0
        fx1[1] <- x1[0]*x1[0] - x1[1] - 5.0
        Fsolve.CopyArray n fx1 fx               // Copy fx1 array values to fx Pointer
        ()                                               // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func3 = Fsolve.FunctionToSolve (funcToSolve3)     // Wrap function so it can be called
    let unknownVariables3 = 2                             // Give number of variables 
    let xGuess3a:double array = [| 0.0; 0.0 |]            // Give a guess value
    // Two separate guess values will be required to obtain both the solutions.
    // This guess will give the solution x = 2, y = -1
    let solveResult3a = Fsolve.Fsolver(func3, unknownVariables3, xGuess3a, Tolerance) // Call solver
    let (soln3a, fx3a, solutionCode3a) = solveResult3a    // Returns solution:
    // soln3a = Array containing solution
    // fx3a = Values of equations at soln3a (should be close to zero within Tolerance)
    // solutionCode3a = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln3a 8                // Prints the solution to '8' decimals
     
    let xGuess3b:double array = [| -2.0; 0.0 |]           // Give another guess value to obtain second solution        
                                                          // This guess will give the solution x = -3, y = 4
    let solveResult3b = Fsolve.Fsolver(func3, unknownVariables3, xGuess3b, Tolerance) 
    let (soln3b, fx3b, solutionCode3b) = solveResult3b
    // soln3b = Array containing solution
    // fx3b = Values of equations at soln3b (should be close to zero within Tolerance)
    // solutionCode3b = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln3b 8                // Prints the solution to '8' decimals



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
        Fsolve.CopyArray n fx1 fx               // Copy fx1 array values to fx Pointer
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
    // e**(−e**(−(x1+x2))) = x2*(1+x1**2)
    // x1*cos(x2) + x2*sin(x1) = 0.5
    // 
    // Solutions:
    // 0.3532    0.6061
   
    // STEP 1: 
    // Define the callback function for this system
    let funcToSolve5 (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature
        let x1 = Fsolve.MakeArray n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx    // Make an array for 'fx' equation/function values from its Pointer                   
        fx1[0] <- exp(-exp(-(x1[0]+x1[1]))) - x1[1]*(1.0+x1[0]**2) // Write equations/functions as f(x) = 0
        fx1[1] <- x1[0]*cos x1[1] + x1[1]*sin x1[0] - 0.5
            
        Fsolve.CopyArray n fx1 fx               // Copy fx1 array values to fx Pointer
        ()                                               // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func5 = Fsolve.FunctionToSolve (funcToSolve5)      // Wrap function so it can be called
    let unknownVariables5 = 2                              // Give number of variables 
    let xGuess5:double array = Array.zeroCreate 2          // Give a guess value
    let solveResult5 = Fsolve.Fsolver(func5, unknownVariables5, xGuess5, Tolerance) // Call solver
    let (soln5, fx5, solutionCode5) = solveResult5        // Returns solution:
    // soln5 = Array containing solution
    // fx5 = Values of equations at soln5 (should be close to zero within Tolerance)
    // solutionCode5 = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln5 4                  // Prints the solution to '7' decimals



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
    // 2x1 − x2  = e**(−x1)
    // −x1 + 2x2 = e**(−x2)
    //
    // Solutions:
    // 0.5671    0.5671
   
    // STEP 1: 
    // Define the callback function for this system
    let funcToSolve6 (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature
        let x1 = Fsolve.MakeArray n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx    // Make an array for 'fx' equation/function values from its Pointer                   
        fx1[0] <- 2.0*x1[0] - x1[1] - exp(-x1[0])        // Write equations/functions as f(x) = 0
        fx1[1] <- -x1[0] + 2.0*x1[1] - exp(-x1[1])
        Fsolve.CopyArray n fx1 fx               // Copy fx1 array values to fx Pointer
        ()                                               // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func6 = Fsolve.FunctionToSolve (funcToSolve6)      // Wrap function so it can be called
    let unknownVariables6 = 2                              // Give number of variables 
    let xGuess6:double array = Array.zeroCreate 2          // Give a guess value
    let solveResult6 = Fsolve.Fsolver(func6, unknownVariables6, xGuess6, Tolerance) // Call solver
    let (soln6, fx6, solutionCode6) = solveResult6        // Returns solution:
    // soln6 = Array containing solution
    // fx6 = Values of equations at soln6 (should be close to zero within Tolerance)
    // solutionCode6 = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln6 4                  // Prints the solution to '7' decimals



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
    // Define the callback function for this system
    let funcToSolve7 (n:int) (x: IntPtr) (fx: IntPtr) =  // This is the function signature
        let x1 = Fsolve.MakeArray n x      // Make an array for 'x' values from its Pointer
        let fx1 = Fsolve.MakeArray n fx    // Make an array for 'fx' equation/function values from its Pointer           
        let a = 1.0
        let b = 1.0
        let c = 1.0
        fx1[0] <- a*(1.0+cos(b))**2 - 4.0*x1[0] * exp(-2.0*c*(a-x1[0])**2) // Write equations/functions as f(x) = 0
        Fsolve.CopyArray n fx1 fx               // Copy fx1 array values to fx Pointer
        ()                                               // Returns equivalent of void in 'C'
    
    // STEP 2:
    // Solve the function
    let func7 = Fsolve.FunctionToSolve (funcToSolve7)      // Wrap function so it can be called
    let unknownVariables7 = 1                              // Give number of variables 
    let xGuess7:double array = Array.zeroCreate 1          // Give a guess value
    let solveResult7 = Fsolve.Fsolver(func7, unknownVariables7, xGuess7, Tolerance) // Call solver
    let (soln7, fx7, solutionCode7) = solveResult7        // Returns solution:
    // soln7 = Array containing solution
    // fx7 = Values of equations at soln6 (should be close to zero within Tolerance)
    // solutionCode7 = String providing information on exit code
    Fsolve.PrintArray "xSolution" soln7 8                  // Prints the solution to '7' decimals

    0
