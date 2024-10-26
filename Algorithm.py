def matrix_multiply(A, B):
    # Check if B is a vector (list) or matrix (list of lists)
    if isinstance(B[0], list):  # B is a matrix
        cols_B = len(B[0])
    else:  # B is a vector
        cols_B = 1
        B = [[b_elem] for b_elem in B]  # Convert to column matrix

    rows_A = len(A)
    cols_A = len(A[0])

    rows_B = len(B)
    cols_B = len(B[0])

    if cols_A != rows_B:
        raise ValueError("Matrix A's columns must match Matrix B's rows.")

    # Initialize the result matrix with zeroes
    result = [[0] * cols_B for _ in range(rows_A)]

    # Perform matrix multiplication
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                result[i][j] += A[i][k] * B[k][j]

    # If the result is a matrix with one column, flatten it to a vector
    if cols_B == 1:
        result = [result[i][0] for i in range(rows_A)]

    return result


def transpose(A):
    if isinstance(A[0], list):
        return [list(i) for i in zip(*A)]
    else:
        # A is a vector
        return [[elem] for elem in A]


def diag(x):
    return [[x[i] if i == j else 0 for j in range(len(x))] for i in range(len(x))]


def identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def matrix_sum(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def matrix_difference(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def inverse_matrix(A):
    # Using Gauss-Jordan elimination
    n = len(A)
    # Create an augmented matrix [A | I]
    augmented = [A[i] + identity(n)[i] for i in range(n)]

    # Forward elimination
    for i in range(n):
        # Make the diagonal contain all ones
        factor = augmented[i][i]
        if factor == 0:
            raise ValueError("Matrix is singular and cannot be inverted.")
        for j in range(2 * n):
            augmented[i][j] /= factor

        # Make the other elements in column i zeros
        for k in range(n):
            if k != i:
                factor = augmented[k][i]
                for j in range(2 * n):
                    augmented[k][j] -= factor * augmented[i][j]

    # Extract the inverse matrix
    inverse = [row[n:] for row in augmented]
    return inverse


def ones(n, k):
    return [[1 for _ in range(k)] for _ in range(n)]


def multiply_by_num(num, matrix):
    return [[num * matrix[i][j] for j in range(len(matrix[0]))] for i in range(len(matrix))]


def max_abs_difference(x_new, x_old):
    return max([abs(x_new[i][0] - x_old[i][0]) for i in range(len(x_new))])


def check_feasibility(A, x, b):
    Ax = matrix_multiply(A, x)
    for i in range(len(Ax)):
        if abs(Ax[i] - b[i][0]) > 1e-6:
            return False
    for xi in x:
        if xi[0] <= 0:
            return False
    return True


c = [1, 1, 0, 0]
c = transpose([c])
A = [[2, 4, 1, 0], [1, 3, 0, -1]]
b = [16, 9]
b = transpose([b])
x_0 = [1, 1, 10, -5]
# x_0 = [-9, 6, 10, 0]
alpha = 0.5
epsilon = 1e-6  # Approximation accuracy

print(c, A, b, transpose([x_0]), alpha)
print(matrix_multiply(A, x_0))

# Ensure x is a column vector
x = [[xi] for xi in x_0]
# Ensure c is a column vector
c = [[ci] for ci in c]

if matrix_multiply(A, x) == b or check_feasibility(A, x, b):
    iteration = 0
    solved = False
    
    while not solved and iteration < 100:
        iteration += 1
        print(f"Iteration {iteration}")
        
        # Step 1: D = diag(x)
        x_vector = [xi[0] for xi in x]
        D = diag(x_vector)
        # Step 2: Compute A_tilde = A * D
        A_tilde = matrix_multiply(A, D)
        # Step 3: Compute c_tilde = D * c
        c_tilde = matrix_multiply(D, c)
        # Step 4: Compute P = I - A_tilde^T * (A_tilde * A_tilde^T)^-1 * A_tilde
        A_tilde_T = transpose(A_tilde)
        A_tilde_A_tilde_T = matrix_multiply(A_tilde, A_tilde_T)
        try:
            Inv_A_tilde_A_tilde_T = inverse_matrix(A_tilde_A_tilde_T)
        except ValueError as e:
            print(f"Error: {e}")
            break
        Middle_term = matrix_multiply(A_tilde_T, Inv_A_tilde_A_tilde_T)
        Middle_term = matrix_multiply(Middle_term, A_tilde)
        P = matrix_difference(identity(len(x_vector)), Middle_term)
        # Step 5: Compute c_p = P * c_tilde
        c_p = matrix_multiply(P, c_tilde)
        # Step 6: Compute v = |min{c_p_i | c_p_i < 0}|
        c_p_flat = [c_p[i][0] for i in range(len(c_p))]
        negative_c_p = [cpi for cpi in c_p_flat if cpi < 0]
        if negative_c_p:
            v = abs(min(negative_c_p))
        else:
            print("The method is not applicable!")
            break
        # Step 7: Compute x_tilde_new = ones(n, 1) + (alpha / v) * c_p
        scaled_cp = multiply_by_num(alpha / v, c_p)
        ones_matrix = ones(len(c_p), 1)
        x_tilde_new = matrix_sum(ones_matrix, scaled_cp)
        # Step 8: Compute x_new = D * x_tilde_new
        x_new = matrix_multiply(D, x_tilde_new)
        # Step 9: Check convergence
        delta = max_abs_difference(x_new, x)
        if delta <= epsilon:
            solved = True
            x_star = x_new
            print("Optimal solution found:")
            print([x_star[i][0] for i in range(len(x_star))])
            # Compute objective function value
            z = sum([c[i][0] * x_star[i][0] for i in range(len(c))])
            print(f"Objective function value: {z}")
            break
        else:
            x = x_new
    else:
        if not solved:
            print("The method did not converge within the maximum number of iterations.")
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
