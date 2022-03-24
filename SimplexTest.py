import sys

M = 1j
optim = False


# Блок инициализации задачи, выбор базисов и добавление искусственных переменных
def ini(m, n, matrix, Z, alg):
    baz = 3
    if alg:
        Cb = [-M] * m
    else:
        Cb = [M] * m
    basis = [-1] * m
    for i in range(0, n):
        flag = True
        count = 0
        for j in range(0, m):
            if matrix[j][i] == 1:
                indexN = i
                indexM = j
                count += 1
                if count > 1:
                    flag = False
                    break
            elif matrix[j][i] == 0:
                continue
            else:
                flag = False
                break
        if flag:
            Cb[indexM] = Z[indexN]
            basis[indexM] = indexN
            baz -= 1

    for i in range(baz):
        if alg:
            Z.append(-M)
        else:
            Z.append(M)
    k = 0
    for j in [i for i in range(m) if basis[i] == -1]:
        basis[j] = len(Z) - baz + k
        k += 1

    return basis, Cb, Z


# Составление новой таблицы с добавленными искусственными базисами
def fir(m, n, l, matrix, basis):
    table = ny(m, l)
    for i in range(0, m):
        for j in range(l):
            if j < n:
                table[i][j] = matrix[i][j]
            if n <= j == basis[i]:
                table[i][j] = 1
    return table


# Расчет дельта строк fj-cj, вместо числа М используется комплексная единица, проверяется оптимальность плана
# Оптимальность проверяется на все отрицательные|положительные значения в строке дельта
def sec(m, n, Z, Cb, basis, table, flag):
    delt = []
    delt_M = []
    count = 0
    for i in range(n):
        delt.append(0)
        delt_M.append(0)
        for j in range(m):
            if Cb[j] == M or Cb[j] == -M:
                delt_M[i] += Cb[j] * table[j][i]
            else:
                delt[i] += Cb[j] * table[j][i]
        if (Z[i] == M or Z[i] == -M) and i in basis:
            delt_M[i] -= Z[i]
        elif Z[i].imag == 0:
            delt[i] -= Z[i]
    for i in range(0, n):
        if delt_M[i].imag < 0:
            count += 1
        elif delt_M[i].imag == 0 and delt[i] <= 0:
            count += 1
    if not any([i for i in delt_M if i.imag < 0]) and not any([i for i in delt if i < 0]):
        print("Unsuccessful")
        return 0, True
    if count == n:
        print("Получен оптимальный план")
        return 0, True
    print("Дельта = ", delt)
    print("Дельта М = ", delt_M)
    return col(n, delt_M, delt, flag)


# Разрешающий столбец матрицы
def col(n, delt_M, delt, flag):
    mainCol = 0
    if flag:  # MAX
        for j in range(n):
            if delt_M[j].imag < 0 and delt_M[j].imag <= delt_M[mainCol].imag:
                mainCol = j
            elif delt[j] <= delt[mainCol] < 0 and delt_M[j].imag == 0 and delt_M[mainCol].imag == 0:
                mainCol = j
    else:  # MIN
        for j in range(n):
            if delt_M[j].imag > 0 and delt_M[j].imag >= delt_M[mainCol].imag:
                mainCol = j
            elif delt[j] >= delt[mainCol] and delt[j] > 0 and delt_M[j].imag == 0 and delt_M[mainCol].imag == 0:
                mainCol = j
    print("Разрешающий столбец = ", mainCol)
    return mainCol, False


# Определение разрешающей строки
def thir(m, A, table, mainCol):
    mainRow = 0
    mintet = sys.maxsize
    tet = []
    for i in range(m):
        if (table[i][mainCol] >= 0 and A[i] >= 0) or (table[i][mainCol] < 0 and A[i] < 0):
            try:
                tet.append(A[i] / table[i][mainCol])
            except ZeroDivisionError:
                tet.append(-1)
        elif (table[i][mainCol] < 0 and A[i] >= 0) or (table[i][mainCol] >= 0 and A[i] < 0):
            tet.append(-1)
        if 0 <= tet[i] < mintet:
            mainRow = i
            mintet = tet[i]
    print("tet = ", tet)
    print("Разрешающая строка = ", mainRow)
    return mainRow


# Новая таблица
def fou(m, l, mainRow, mainCol, table, A, Cb, Z, basis):
    new_table = ny(m, l)
    newA = []
    for i in range(m):
        if i != mainRow:
            for j in range(l):
                if j != mainCol:
                    new_table[i][j] = ((table[i][j] * table[mainRow][mainCol] - table[i][mainCol] * table[mainRow][j]) /
                                       table[mainRow][mainCol])
                else:
                    new_table[i][mainCol] = 0
            newA.append((A[i] * table[mainRow][mainCol] - A[mainRow] * table[i][mainCol]) / table[mainRow][mainCol])
        else:
            newA.append((A[mainRow] / table[mainRow][mainCol]))
    for j in range(l):
        new_table[mainRow][j] = (
                table[mainRow][j] / table[mainRow][mainCol])
    if Cb[mainRow].imag != 0:
        for i in range(0, m):
            new_table[i][basis[mainRow]] = 0
    Cb[mainRow] = Z[mainCol]
    basis[mainRow] = mainCol
    print("\nНовая матрица")
    print(new_table)
    print("Базисные переменные = ", basis)
    return new_table, newA, Cb, basis


# Заполнение плана
def fif(m, basis, A, l):
    result = []
    for j in range(l):
        result.append(0)
    for i in range(m):
        result[basis[i]] = A[i]
    return result


# Проверка на неограниченность
def StopCriterion(m, table, main):
    count = 0
    for i in range(m):
        if table[i][main] <= 0:
            count += 1
    if count == m:
        print("\nUnbound")
        return False
    return True


# Расчет значения целевой функции
def F(m, Cb, A):
    z = 0
    for i in range(m):
        z += Cb[i] * A[i]
    return z


# Создание заполненной нулями таблицы
def ny(m, n):
    new_table = [0] * m
    for i in range(m):
        new_table[i] = [0] * n
    return new_table


# Последовательная работа алгоритма вычисления
def pack(source, A, Z, isMax):
    m = len(source)
    n = len(source[0])

    basis, Cb, Z = ini(m, n, source, Z, isMax)

    l = len(Z)

    table = fir(m, n, l, source, basis)

    print("Начальная таблица = ", table)
    print("Базисные переменные = ", basis)

    mainCol, optim = sec(m, l, Z, Cb, basis, table, isMax)

    while not optim and StopCriterion(m, table, mainCol):
        mainRow = thir(m, A, table, mainCol)
        table, A, Cb, basis = fou(m, l, mainRow, mainCol, table, A, Cb, Z, basis)
        mainCol, optim = sec(m, l, Z, Cb, basis, table, isMax)

    result = fif(m, basis, A, l)
    print("План = ", result)

    f = F(m, Cb, A)
    print("Значение целевой функции = ", f)
