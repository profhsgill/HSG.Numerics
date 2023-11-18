Imports HSG.Numerics
Module FsolveTester

    Sub Main()

        Dim Tolerance As Double = 0.0000000001 ' Tolerance for all test cases

        ' Solve function 1
        ' ----------------
        Dim func1 As New Fsolve.FunctionToSolve(AddressOf funcToSolve1) ' Wrap function so it can be called
        Dim unknownVariables1 As Integer = 1                            ' Give number of variables 
        Dim xGuess1() = {0.0}                                           ' Give a guess value
        Dim solveResult1 = Fsolve.Fsolver(func1, unknownVariables1, xGuess1, Tolerance) ' Call solver
        ' Returns solution:
        ' soln1 = Array containing solution
        ' fx1 = Values of equations at soln1 (should be close to zero within Tolerance)
        ' solutionCode1 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult1.Item1, 8)


        ' Solve function 2
        ' ----------------
        Dim func2 As New Fsolve.FunctionToSolve(AddressOf funcToSolve2) ' Wrap function so it can be called
        Dim unknownVariables2 As Integer = 3                            ' Give number of variables 
        Dim xGuess2() = {0.0, 0.0, 0.0}                                 ' Give a guess value
        Dim solveResult2 = Fsolve.Fsolver(func2, unknownVariables2, xGuess2, Tolerance) ' Call solver
        ' Returns solution:
        ' soln2 = Array containing solution
        ' fx2 = Values of equations at soln2 (should be close to zero within Tolerance)
        ' solutionCode2 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult2.Item1, 8)



        ' Solve function 3
        ' ----------------
        Dim func3 As New Fsolve.FunctionToSolve(AddressOf funcToSolve3) ' Wrap function so it can be called
        Dim unknownVariables3 As Integer = 2                            ' Give number of variables 
        Dim xGuess3a() = {0.0, 0.0}                                     ' Give a guess value
        ' Two separate guess values will be required to obtain both the solutions.
        ' This guess will give the solution x = 2, y = -1
        Dim solveResult3a = Fsolve.Fsolver(func3, unknownVariables3, xGuess3a, Tolerance) ' Call solver
        ' Returns solution:
        ' soln3a = Array containing solution
        ' fx3a = Values of equations at soln3a (should be close to zero within Tolerance)
        ' solutionCode3a = String providing information on exit code
        Fsolve.PrintArray("x", solveResult3a.Item1, 8)

        ' Give 2nd guess
        '----------------
        Dim xGuess3b() = {-2.0, 0.0}  ' Give a guess value
        ' Two separate guess values will be required to obtain both the solutions.
        ' This guess will give the solution x = -3, y = 4
        Dim solveResult3b = Fsolve.Fsolver(func3, unknownVariables3, xGuess3b, Tolerance) ' Call solver
        ' Returns solution:
        ' soln3b = Array containing solution
        ' fx3b = Values of equations at soln3b (should be close to zero within Tolerance)
        ' solutionCode3a = String providing information on exit code
        Fsolve.PrintArray("x", solveResult3b.Item1, 8)


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
    ' ....... TEST 1 .............................
    ' ....... Linear System ......................
    ' ....... Unknown variables = 1 ..............
    ' ============================================
    '
    ' Equations:
    ' 1) x = 1
    '
    ' Solution:
    ' x = 1
    Sub funcToSolve1(n As Integer, x As IntPtr, fx As IntPtr)         ' This is the function signature
        Dim x1() As Double = Fsolve.ExtractArrayFromPointer(n, x)     ' Make an array for 'x' values from its Pointer
        Dim fx1() As Double = Fsolve.ExtractArrayFromPointer(n, fx)   ' Make an array for 'fx' equation/function values from its Pointer
        fx1(0) = x1(0) - 1.0                                          ' Write equations/functions as f(x) = 0, here it becomes x - 1 = 0
        Fsolve.CopyArrayToPointer(n, fx1, fx)                         ' Copy fx array values to fx Pointer
    End Sub


    ' ============================================
    ' ....... TEST 2 .............................
    ' ....... Linear System ......................
    ' ....... Unknown variables = 3 ..............
    ' ============================================
    '
    ' Equations:
    ' -x + 3y + 7z = 0
    ' 2x - 2y - z = 0
    ' x + y + z = 1
    '
    ' Solution:
    ' x = 0.55, y = 0.65, and z = -0.2
    Sub funcToSolve2(n As Integer, x As IntPtr, fx As IntPtr)         ' This is the function signature
        Dim x1() As Double = Fsolve.ExtractArrayFromPointer(n, x)     ' Make an array for 'x' values from its Pointer
        Dim fx1() As Double = Fsolve.ExtractArrayFromPointer(n, fx)   ' Make an array for 'fx' equation/function values from its Pointer
        fx1(0) = -1.0 * x1(0) + 3.0 * x1(1) + 7.0 * x1(2)             ' Write equations/functions as f(x) = 0
        fx1(1) = 2.0 * x1(0) - 2.0 * x1(1) - 1.0 * x1(2)
        fx1(2) = x1(0) + x1(1) + x1(2) - 1.0
        Fsolve.CopyArrayToPointer(n, fx1, fx)                         ' Copy fx array values to fx Pointer
    End Sub


    ' ============================================
    ' ....... TEST 3 .............................
    ' ....... Non Linear System ..................
    ' ....... Unknown variables = 2 ..............
    ' ============================================
    '
    ' Equations:
    ' x + y = 1
    ' y = x**2 - 5
    '
    ' Solutions:
    ' (1) x = -3, y = 4;      (2) x = 2, y = -1
    Sub funcToSolve3(n As Integer, x As IntPtr, fx As IntPtr)         ' This is the function signature
        Dim x1() As Double = Fsolve.ExtractArrayFromPointer(n, x)     ' Make an array for 'x' values from its Pointer
        Dim fx1() As Double = Fsolve.ExtractArrayFromPointer(n, fx)   ' Make an array for 'fx' equation/function values from its Pointer
        fx1(0) = x1(0) + x1(1) - 1.0                                  ' Write equations/functions as f(x) = 0
        fx1(1) = x1(0) * x1(0) - x1(1) - 5.0
        Fsolve.CopyArrayToPointer(n, fx1, fx)                         ' Copy fx array values to fx Pointer
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