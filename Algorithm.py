def matrix_multiply(A, B):
    # Get the number of rows and columns of matrix A and matrix B
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    # Fix: Check if B is a list of lists for proper matrix operation (column vector handling)
    cols_B = len(B[0]) if isinstance(B[0], list) else 1

    # Check if matrices can be multiplied (columns of A should equal rows of B)
    if cols_A != rows_B:
        raise ValueError("Matrix A's columns must match Matrix B's rows.")

    # Initialize the result matrix with zeroes
    result = [[0] * cols_B for _ in range(rows_A)]

    # Perform matrix multiplication
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                # Fix: Proper multiplication check for column vector multiplication
                if not isinstance(B[0], list):  # Case for matrix * column vector
                    result[i][j] += A[i][k] * B[k]
                else:  # Case for matrix * matrix
                    result[i][j] += A[i][k] * B[k][j]

    return result


def transpose(A):
    # Fix: Ensure correct row and column swapping in transpose
    num_rows = len(A)
    num_cols = len(A[0]) if isinstance(A[0], list) else 1
    result = [[0] * num_rows for _ in range(num_cols)]

    if isinstance(A[0], list):
        for i in range(num_rows):
            for j in range(num_cols):
                result[j][i] = A[i][j]
    else:
        # Fix: If it's a vector (1D array), transpose directly to column
        for i in range(num_cols):
            result[i][0] = A[i]

    return result


def diag(x):
    answer = [[0] * len(x) for _ in range(len(x))]
    for i in range(len(x)):
        answer[i][i] = x[i]
    return answer


def identity(n):
    answer = [[0] * n for _ in range(n)]
    for i in range(len(answer)):
        answer[i][i] = 1
    return answer


def matrix_sum(A, B):
    # Check if the dimensions match for addition
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError("Matrices must have the same dimensions for addition.")

    # Initialize result matrix
    result = [[0 for _ in range(len(A[0]))] for _ in range(len(A))]

    # Perform element-wise addition
    for i in range(len(A)):
        for j in range(len(A[0])):
            result[i][j] = A[i][j] + B[i][j]

    return result


def matrix_difference(A, B):
    # Check if the dimensions match for subtraction
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        raise ValueError("Matrices must have the same dimensions for subtraction.")

    # Initialize result matrix
    result = [[0 for _ in range(len(A[0]))] for _ in range(len(A))]

    # Perform element-wise subtraction
    for i in range(len(A)):
        for j in range(len(A[0])):
            result[i][j] = A[i][j] - B[i][j]

    return result


def inverse_matrix(A):
    n = len(A)

    # Create the identity matrix of the same size as A
    I = [[float(i == j) for i in range(n)] for j in range(n)]

    # Augment the matrix A with the identity matrix I
    augmented_matrix = [A[i] + I[i] for i in range(n)]

    # Perform Gaussian elimination
    for i in range(n):
        # Check if the diagonal element is zero (to avoid division by zero)
        if augmented_matrix[i][i] == 0:
            raise ValueError("Matrix is singular and cannot be inverted.")

        # Normalize the current row by dividing by the diagonal element
        diag_element = augmented_matrix[i][i]
        for j in range(2 * n):
            augmented_matrix[i][j] /= diag_element

        # Eliminate the current column for all other rows
        for k in range(n):
            if k != i:
                factor = augmented_matrix[k][i]
                for j in range(2 * n):
                    augmented_matrix[k][j] -= factor * augmented_matrix[i][j]

    # Extract the inverse matrix (the right half of the augmented matrix)
    inverse = [row[n:] for row in augmented_matrix]

    return inverse


def ones(n, k):
    answer = [[0 for _ in range(k)] for _ in range(n)]
    for i in range(n):
        for j in range(k):
            answer[i][j] = 1
    return answer


def multiply_by_num(num, matrix):
    answer = [[0 for _ in range(len(matrix[0]))] for _ in range(len(matrix))]
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            answer[i][j] = num * matrix[i][j]
    # print(answer)
    return answer


def check_condition(A, x_0, b):
    matrix = matrix_difference(matrix_multiply(A, x_0), b)
    print(matrix)
    maximal = -1
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if abs(matrix[i][j]) > maximal:
                maximal = abs(matrix[i][j])

    if maximal <= 0.000001:
        return True
    else:
        return False


c = [1, 1, 0, 0]
c = transpose([c])  # Fix: Transpose properly into a column vector
A = [[2, 4, 1, 0], [1, 3, 0, -1]]
b = [16, 9]
b = transpose([b])
# x_0 = [1, 1, 10, -5]
x_0 = [-9, 6, 10, 0]
alpha = 0.5
print(c, A, b, transpose([x_0]), alpha)
print(matrix_multiply(A, x_0))
if matrix_multiply(A, x_0) == b or check_condition(A, x_0, b): # or max(abs(x) for x in matrix_difference(matrix_multiply(A, transpose([x_0])), b)) <= 1e-6)
    x = x_0
    solved = False
    while not solved:
        D = diag(x)
        print("D:", D)
        A_telda = matrix_multiply(A, D)
        print("A_telda:", A_telda)
        c_telda = matrix_multiply(D, c)
        print("c_telda:", c_telda)
        P = matrix_difference(identity(len(x)), matrix_multiply(matrix_multiply(transpose(A_telda), inverse_matrix(matrix_multiply(A_telda, transpose(A_telda)))), A_telda))
        print("P1:", P)
        # P = matrix_difference(identity(len(x)), P)
        # print("P2:", P)
        c_p = matrix_multiply(P, c_telda)
        print("c_p:", c_p)
        c_p_flat = [c_p[i][0] for i in range(len(c_p))]
        print("c_p_flat:", c_p_flat)
        negative_numbers = [x for x in c_p_flat if x < 0]
        print("negative_numbers:", negative_numbers)

        # Find the minimal negative number and convert it to positive
        if len(negative_numbers) > 0:
            v = abs(min(negative_numbers))
        else:
            v = None

        if v is not None:
            print("v:", v)
            x_telda = matrix_sum(ones(len(c_p), 1), multiply_by_num(alpha / v,  c_p))
            print("x_telda:", x_telda)
            x_star = matrix_multiply(D, x_telda)
            print("x_star:", x_star)
            # if max(abs(x) for x in matrix_difference(x_telda, x_star)) < 1e-6:
            print(matrix_multiply(identity(len(x_telda)), x_telda), x_telda)
            if check_condition(identity(len(x)), x, x_star):
                solved = True
                print(x_star)
            x = [x_star[i][0] for i in range(len(x_star))]
        else:
            print("Method is not applicable!")
else:
    print("Not correct initial solution")

# if (A * x_0 == b) or (max(abs(A * x_0 - b)) <= 0.00001
# start solve
#   x = x_0
#   solved = False
#   while solved == False:
#       D = diag(x)
#       A_telda = A * D
#       c_telda = D * c
#       P = eye(n) - A_telda' * (A_telda * A_telda')^(-1) * A_telda
#       P = eye(length(P)) - P
#       c_p = P * c-telda
#       v = abs(min(c_p(c_p < 0))) - choose minimal negative number from c_p and convert it to positive
#       x_telda = ones(length(c_p), 1) + (alpha_/v) * c_p
#       x_start = D * x_telda
#       if max(abs(x_telda - x_star)) <= 0.00001:
#           solved = True
#           print(x_star)
#       x = x_star
# else
#   print("Not correct initial solution")
