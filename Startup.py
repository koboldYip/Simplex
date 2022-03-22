import SimplexTest
from scipy.optimize import linprog

obj = [4, -2, 3, 0]
lhs_uneq = [[-3, -2, -1, 0]]
rhs_uneq = [0]
lhs_eq = [[-2, -1, -3, 1], [-4, -4, 3, 1]]
rhs_eq = [4, 6]

print(linprog(c=obj, A_ub=lhs_uneq, b_ub=rhs_uneq,
              A_eq=lhs_eq, b_eq=rhs_eq,
              method="revised simplex"))

print()

srce = [[-2, -1, -3, 1, 0], [-4, -4, 3, 1, 0], [-3, -2, -1, 0, 1]]  # Левые части уравнений
A = [4, 6, 0]  # Правые части уравнений
alg = True  # max = True; min = False
Z = [-4, 2, -3, 0, 0]

SimplexTest.pack(srce, A, Z, alg)
