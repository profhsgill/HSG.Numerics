Imports HSG.Numerics
Module FsolveTester
    Sub Main()

        Dim Tolerance As Double = 0.0000000001 ' Tolerance for all test cases

        ' Solve function 1
        ' ----------------
        Dim unknownVariables1 As Integer = 1                            ' Give number of variables 
        Dim xGuess1() = {0.0}                                           ' Give a guess value
        Dim solveResult1 = Fsolve.Fsolver(AddressOf funcToSolve1, unknownVariables1, xGuess1, Tolerance) ' Call solver
        ' Returns solution:
        ' soln1 = Array containing solution
        ' fx = Values of equations at soln1 (should be close to zero within Tolerance)
        ' solutionCode1 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult1.Item1, 8)



        ' Solve function 2
        ' ----------------
        Dim unknownVariables2 As Integer = 3                            ' Give number of variables 
        Dim xGuess2() = {0.0, 0.0, 0.0}                                 ' Give a guess value
        Dim solveResult2 = Fsolve.Fsolver(AddressOf funcToSolve2, unknownVariables2, xGuess2, Tolerance) ' Call solver
        ' Returns solution:
        ' soln2 = Array containing solution
        ' fx2 = Values of equations at soln2 (should be close to zero within Tolerance)
        ' solutionCode2 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult2.Item1, 8)



        ' Solve function 3
        ' ----------------
        Dim unknownVariables3 As Integer = 2                            ' Give number of variables 
        Dim xGuess3a() = {0.0, 0.0}                                     ' Give a guess value
        ' Two separate guess values will be required to obtain both the solutions.
        ' This guess will give the solution x = 2, y = -1
        Dim solveResult3a = Fsolve.Fsolver(AddressOf funcToSolve3, unknownVariables3, xGuess3a, Tolerance) ' Call solver
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
        Dim solveResult3b = Fsolve.Fsolver(AddressOf funcToSolve3, unknownVariables3, xGuess3b, Tolerance) ' Call solver
        ' Returns solution:
        ' soln3b = Array containing solution
        ' fx3b = Values of equations at soln3b (should be close to zero within Tolerance)
        ' solutionCode3a = String providing information on exit code
        Fsolve.PrintArray("x", solveResult3b.Item1, 8)


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


        ' Solve function 5
        ' ----------------
        Dim unknownVariables5 As Integer = 2                            ' Give number of variables 
        Dim xGuess5() = {0.0, 0.0}   ' Give a guess value
        Dim solveResult5 = Fsolve.Fsolver(AddressOf funcToSolve5, unknownVariables5, xGuess5, Tolerance) ' Call solver
        ' Returns solution:
        ' soln5 = Array containing solution
        ' fx5 = Values of equations at soln5 (should be close to zero within Tolerance)
        ' solutionCode4 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult5.Item1, 4)


        ' Solve function 6
        ' ----------------
        Dim unknownVariables6 As Integer = 2                            ' Give number of variables 
        Dim xGuess6() = {0.0, 0.0}   ' Give a guess value
        Dim solveResult6 = Fsolve.Fsolver(AddressOf funcToSolve6, unknownVariables6, xGuess6, Tolerance) ' Call solver
        ' Returns solution:
        ' soln6 = Array containing solution
        ' fx6 = Values of equations at soln6 (should be close to zero within Tolerance)
        ' solutionCode6 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult6.Item1, 4)


        ' Solve function 7
        ' ----------------
        Dim unknownVariables7 As Integer = 1                            ' Give number of variables 
        Dim xGuess7() = {0.0}   ' Give a guess value
        Dim solveResult7 = Fsolve.Fsolver(AddressOf funcToSolve7, unknownVariables7, xGuess7, Tolerance) ' Call solver
        ' Returns solution:
        ' soln7 = Array containing solution
        ' fx7 = Values of equations at soln7 (should be close to zero within Tolerance)
        ' solutionCode7 = String providing information on exit code
        Fsolve.PrintArray("x", solveResult7.Item1, 8)

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
    Function funcToSolve1(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        fx(0) = x(0) - 1.0 ' Write equations/functions as f(x) = 0, here it becomes x - 1 = 0
        Return fx
    End Function


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
    Function funcToSolve2(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        fx(0) = -1.0 * x(0) + 3.0 * x(1) + 7.0 * x(2)             ' Write equations/functions as f(x) = 0
        fx(1) = 2.0 * x(0) - 2.0 * x(1) - 1.0 * x(2)
        fx(2) = x(0) + x(1) + x(2) - 1.0
        Return fx
    End Function


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
    Function funcToSolve3(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        fx(0) = x(0) + x(1) - 1.0                                  ' Write equations/functions as f(x) = 0
        fx(1) = x(0) * x(0) - x(1) - 5.0
        Return fx
    End Function


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
    ' (3 � 2*x(0)) * x(0)                     -2*x(1)                             = -1
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



    ' ============================================
    ' ....... TEST 5 .............................
    ' ....... Non Linear System ..................
    ' ....... Unknown variables = 2 ..............
    ' ============================================
    '
    ' This is an example from https://www.mathworks.com/help/optim/ug/fsolve.html
    '
    ' Equations:
    '
    ' e^(-e^(-(x+x2))) = x2*(1+x^2)
    ' x*cos(x2) + x2*sin(x) = 0.5
    ' 
    ' Solutions:
    ' 0.3532    0.6061
    Function funcToSolve5(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        fx(0) = Math.Exp(-Math.Exp(-(x(0) + x(1)))) - x(1) * (1.0 + x(0) ^ 2) ' Write equations/functions as f(x) = 0
        fx(1) = x(0) * Math.Cos(x(1)) + x(1) * Math.Sin(x(0)) - 0.5
        Return fx
    End Function



    ' ============================================
    ' ....... TEST 6 .............................
    ' ....... Non Linear System ..................
    ' ....... Unknown variables = 2 ..............
    ' ============================================
    '
    ' This is an example from https://www.mathworks.com/help/optim/ug/fsolve.html
    '
    ' Equations:
    '
    ' 2x1 - x2  = e^(-x)
    ' -x + 2x2 = e^(-x2)
    '
    ' Solutions:
    ' 0.5671    0.5671
    Function funcToSolve6(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        fx(0) = 2.0 * x(0) - x(1) - Math.Exp(-x(0))               ' Write equations/functions as f(x) = 0
        fx(1) = -x(0) + 2.0 * x(1) - Math.Exp(-x(1))
        Return fx
    End Function



    ' ============================================
    ' ....... TEST 7 .............................
    ' ....... Non Linear System ..................
    ' ....... Unknown variables = 1 ..............
    ' ============================================
    '
    ' This is an example from https://www.mathworks.com/help/optim/ug/fsolve.html
    '
    ' Equations:
    '
    ' a(1+cos(b))^2 = 4x * exp(-2c*(a-x)^2)
    '
    'Depends on a, b, c
    ' 0.70548923 for a = b = c = 1
    Function funcToSolve7(x() As Double) As Double() ' This is the function signature
        Dim n As Integer = x.Length ' Find number of elements in 'x' array 
        Dim fx(n - 1) As Double ' VB uses index to declare array, so use 'n-1' since index starts at '0'
        Dim a As Double = 1.0
        Dim b As Double = 1.0
        Dim c As Double = 1.0
        fx(0) = a * (1.0 + Math.Cos(b)) ^ 2 - 4.0 * x(0) * Math.Exp(-2.0 * c * (a - x(0)) ^ 2)   ' Write equations/functions as f(x) = 0
        Return fx
    End Function


End Module