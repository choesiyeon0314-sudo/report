#행렬 출력
def matrixout(mx):
    n = len(mx)
    print("┌" + "       " * n + "┐")
    for i in range(n):
        print("|", end=" ")
        for j in range(n):
            print("%.3f" % mx[i][j], end=" ")
        print("|")
    print("└" + "       " * n + "┘")

#행렬 입력받기
def input_matrix(n):
    matrix = []
    for i in range(n):
        row = input(f"{i+1}행 (공백으로 구분합니다.): ").strip().split()
        if len(row) != n:
            print(f"입력 오류: 정확히 {n}개의 실수를 입력해야 합니다.")
            return None
        matrix.append([float(x) for x in row])
    return matrix

#전치행렬
def transpose(m):
    return [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]

#가우스조던과 행렬식 방식으로 구한 역행렬이 같은지 비교
def compare_matrices(matrix1, matrix2):
    n = len(matrix1)
    for i in range(n):
        for j in range(n):
            if abs(matrix1[i][j] - matrix2[i][j]) > 1e-12:
                return False
    return True

# 행렬식기반
#소행렬식구하기
def minor_matrix(m, i, j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

#determinant구하기
def determinant_matrix(m):
    n = len(m)
    if n == 1:
        return m[0][0]
    if n == 2:
        return m[0][0]*m[1][1] - m[0][1]*m[1][0]

    determinant = 0
    for c in range(n):
        determinant += ((-1)**c) * m[0][c] * determinant_matrix(minor_matrix(m, 0, c))
    return determinant

#역행렬 - 행렬식
def inverse_matrix(m):
    n = len(m)
    determinant = determinant_matrix(m)
    if abs(determinant) < 1e-12:
        return None  # 역행렬 없음

    if n == 1:
        return [[1.0 / m[0][0]]]
    if n == 2:
        return [
            [ m[1][1]/determinant, -1*m[0][1]/determinant],
            [-1*m[1][0]/determinant,  m[0][0]/determinant]
        ]

    cofactor = []
    for r in range(n):
        row = []
        for c in range(n):
            minor = minor_matrix(m, r, c)
            row.append(((-1)**(r + c)) * determinant_matrix(minor))
        cofactor.append(row)
    adjust = transpose(cofactor)
    for i in range(n):
        for j in range(n):
            adjust[i][j] /= determinant
    return adjust

# ---------- 가우스 조던 ----------
def unit_matrix(n):
    I = [[0.0]*n for _ in range(n)]
    for i in range(n):
        I[i][i] = 1.0
    return I

def temp_matrix_row(m):
    return [row[:] for row in m]

def inverse_gauss_jordan(m):
    n = len(m)
    pri_matrix = temp_matrix_row(m)
    unit_M = unit_matrix(n)

    for col in range(n):
        # 피벗 행 선택
        max_value = 0
        select_row = col
        for r in range(col, n):
            if abs(pri_matrix[r][col]) > max_value:
                max_value = abs(pri_matrix[r][col])
                select_row = r

        # 피벗이 0이면 역행렬 없음
        if abs(pri_matrix[select_row][col]) < 1e-12:
            return None

        # 행 교환
        if select_row != col:
            pri_matrix[col], pri_matrix[select_row] = pri_matrix[select_row], pri_matrix[col]
            unit_M[col], unit_M[select_row] = unit_M[select_row], unit_M[col]

        # 주대각 1로 만들기
        pivot = pri_matrix[col][col]
        for j in range(n):
            pri_matrix[col][j] /= pivot
            unit_M[col][j] /= pivot

        # 주대각 외 성분을 0으로 만들기
        for r in range(n):
            if r == col:
                continue
            factor = pri_matrix[r][col]
            if abs(factor) > 0:
                for j in range(n):
                    pri_matrix[r][j] -= factor * pri_matrix[col][j]
                    unit_M[r][j] -= factor * unit_M[col][j]
    return unit_M


#내가 만든 함수

def show(m, I):
  n = len(m)
  print("-" * (n * 15))
  print("┌" + "       " * 3*n + "┐")
  for i in range(n):
      m_to_i = " ".join(f"{m[i][j]:8.3f}" for j in range(n))
      i_to_inverse = " ".join(f"{I[i][j]:8.3f}" for j in range(n))
      print(f"| {m_to_i} | {i_to_inverse} |")
  print("└" + "       " * 3*n + "┘")


def process_gauss_jordan(m):
    n = len(m)
    m = [row[:] for row in m]
    unit_m = [[float(i == j) for j in range(n)] for i in range(n)]

    print("초기 행렬 [A | I]:")
    show(m, unit_m)

    for col in range(n):
        print(f"\n=== {col+1}번째 열 처리 ===")

        select_row = col
        max_val = 0
        for r in range(col, n):
            if abs(m[r][col]) > max_val:
                max_val = abs(m[r][col])
                select_row = r

        if select_row != col:
            print(f"→ 행 {col+1} ↔ 행 {select_row+1} 교환")
            m[col], m[select_row] = m[select_row], m[col]
            unit_m[col], unit_m[select_row] = unit_m[select_row], unit_m[col]
            show(m, unit_m)

        pivot = m[col][col]
        print(f"→ 피벗 = {pivot:.3f} → 행 {col+1} 전체를 {pivot:.3f}로 나눔")
        for j in range(n):
            m[col][j] /= pivot
            unit_m[col][j] /= pivot
        show(m, unit_m)

        for r in range(n):
            if r == col:
                continue
            factor = m[r][col]
            if abs(factor) > 0:
                print(f"→ 행 {r+1} = 행 {r+1} - ({factor:.3f}) × 행 {col+1}")
                for j in range(n):
                    m[r][j] -= factor * m[col][j]
                    unit_m[r][j] -= factor * unit_m[col][j]
                show(m, unit_m)

    print(" [I | A⁻¹]:")
    show(m, unit_m)



def main():
    n = input("정방행렬의 차수를 입력하세요(n은 1이상의 정수): ").strip()
    if not n.isdigit():
        print("입력 오류: 정수를 입력하세요.")
        return
    n = int(n)

    print(f"\n{n}×{n} 정방행렬 A 를 입력하세요.")
    matrix = input_matrix(n)
    if matrix is None:
        return

    # 행렬식 방식
    inv_adj = inverse_matrix(matrix)
    if inv_adj is None:
        print("\n행렬식이 0이기 때문에 계산 불가합니다.")
    else:
        print("\n행렬식으로 구한 역행렬:")
        matrixout(inv_adj)

    # 가우스-조던 방식
    inv_gj = inverse_gauss_jordan(matrix)
    if inv_gj is None:
        print("\n가우스-조던 계산 불가.")
    else:
        print("가우스-조던 소거법으로 구한 역행렬:")
        matrixout(inv_gj)

    # 결과 비교
    if inv_adj is not None and inv_gj is not None:
        compare_result = compare_matrices(inv_adj, inv_gj)
        if compare_result is True:
            print("\n두 결과가 동일합니다.")
        else:
            print("\n두 결과가 서로 다릅니다.")
    else:
        print("\n행렬식, 가우스-조던으로 역행렬을 구하지 못했습니다.")

    determinant=determinant_matrix(matrix)
    if determinant>0:
      process_gauss_jordan(matrix)

if __name__ == "__main__":
    main()
