#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
const double limit = 0.000001;
class Matrix{
protected:
    int rows, columns;
    double **data{};
    bool change = true;
public:
    Matrix() {
        rows = 0;
        columns = 0;
    }
    Matrix(int row, int column) {
        rows = row;
        columns = column;
        data = new double* [rows];
        for (int i = 0; i < rows; i++) {
            data[i] = new double[columns];
        }
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                data[i][j] = 0.0;
            }
        }
    }
    void printMatrix(){
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                cout.precision(4);
                cout << fixed << data[i][j] << " ";
            }
            cout << endl;
        }
    }
    void setElement(int row, int column, double elem){
        this->data[row][column] = elem;
    }
    int getRows(){
        return rows;
    }
    int getColumns(){
        return columns;
    }
    double getElement(int row, int column){
        return data[row][column];
    }
    Matrix operator+(const Matrix &);
    Matrix operator-(const Matrix &);
    Matrix& operator=(Matrix otherMatrix);
    Matrix operator*(const Matrix &);
    int findMaxColumn(int j){
        double max = 0;
        int index = 0;
        for (int i = j; i < rows; i++){
            if (abs(data[i][j]) > max){
                max = abs(data[i][j]);
                index = i;
            }
        }
        return index;
    }
    static Matrix transpose(const Matrix &);
    friend istream& operator >>(istream&in, Matrix &matrix);
    friend ostream& operator <<(ostream&out, Matrix &matrix);
};
istream& operator>>(istream&in, Matrix &matrix){
    for(int i = 0; i < matrix.rows; i++){
        for (int j = 0; j < matrix.columns; j++){
            in >> matrix.data[i][j];
        }
    }
    return in;
}
ostream& operator<<(ostream&out, Matrix &matrix) {
    for(int i = 0; i < matrix.rows; i++){
        for (int j = 0; j < matrix.columns; j++){
            out.precision(4);
            out << fixed << matrix.data[i][j] <<" ";
        }
        out << endl;
    }
    return out;
}

Matrix Matrix::operator+(const Matrix &otherMatrix) {
    if (otherMatrix.rows != rows || otherMatrix.columns != columns) {
        this->change = false;
        cout << "Error: the dimensional problem occurred\n";
        return *this;
    } else {
        Matrix temp(this->rows, this->columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp.setElement(i, j, this->data[i][j] + otherMatrix.data[i][j]);
            }
        }
        this->change = true;
        return temp;
    }
}
Matrix Matrix::operator-(const Matrix &otherMatrix) {
    if (otherMatrix.rows != rows || otherMatrix.columns != columns) {
        cout << "Error: the dimensional problem occurred\n";
        this->change = false;
        return *this;
    } else {
        Matrix temp(this->rows, this ->columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp.setElement(i, j, this->data[i][j] - otherMatrix.data[i][j]);
            }
        }
        this->change = true;
        return temp;
    }

}
Matrix& Matrix::operator=(Matrix otherMatrix) {
    if (this->rows == otherMatrix.rows && this->columns == otherMatrix.columns){
        for (int i = 0; i < this->rows; i++){
            for (int j = 0; j < this->columns; j++){
                if (fabs(otherMatrix.data[i][j]) < limit){
                    this->data[i][j] = 0.0;
                } else {
                    this->data[i][j] = otherMatrix.data[i][j];
                }
            }
        }
    }
    return *this;
}
Matrix Matrix::operator*(const Matrix &otherMatrix){
    if (this->columns == otherMatrix.rows){
        Matrix temp(this->rows, otherMatrix.columns);
        for (int i = 0; i < this->rows; i++){
            for (int j = 0; j < otherMatrix.columns; j++){
                for (int k = 0; k < this->columns; k++){
                    temp.setElement(i,j, temp.getElement(i,j) + this->data[i][k] * otherMatrix.data[k][j]);
                }
            }
        }
        this->change = true;
        return temp;
    } else{
        cout << "Error: the dimensional problem occurred\n";
        this->change = false;
        return *this;
    }
}

Matrix Matrix::transpose(const Matrix &otherMatrix) {
    Matrix temp(otherMatrix.columns, otherMatrix.rows);
    for(int i = 0; i < otherMatrix.rows; i++){
        for (int j = 0; j < otherMatrix.columns; j++){
            temp.setElement(j, i,  otherMatrix.data[i][j]);
        }
    }
    return temp;
}
class squaredMatrix: public Matrix{
public:
    squaredMatrix(){
        rows = columns = 0;
    }
    explicit squaredMatrix(int n){
        rows = columns = n;
        data = new double* [rows];
        for (int i = 0; i < n; i++){
            data[i] = new double[columns];
        }
        for (int i = 0; i < n; i ++){
            for (int j = 0; j < n; j++){
                data[i][j] = 0;
            }
        }
    }
    squaredMatrix operator+(auto square){
        Matrix* first = this;
        Matrix* second  = &square;
        Matrix multiplication = *first + *second;
        auto *ans = (squaredMatrix*) &multiplication;
        return *ans;
    }
    squaredMatrix operator*(auto square){
        Matrix* first = this;
        Matrix* second  = &square;
        Matrix multiplication = *first * *second;
        auto *ans = (squaredMatrix*) &multiplication;
        return *ans;
    }
    squaredMatrix operator-(auto square){
        Matrix* first = this;
        Matrix* second  = &square;
        Matrix minus = *first - *second;
        auto *ans = (squaredMatrix*) &minus;
        return *ans;

    }
    squaredMatrix operator=(auto square){
        Matrix* first = this;
        Matrix* second  = &square;
        *first = *second;
        auto *ans = (squaredMatrix*) &first;
        return *ans;
    }
    static squaredMatrix transpose(auto square){
        Matrix* matrix = &square;
        Matrix trans = trans.transpose(*matrix);
        auto *ans = (squaredMatrix*) &trans;
        return *ans;
    }
};
class IdentityMatrix: public squaredMatrix{
public:
    IdentityMatrix(){
        rows = columns = 0;
    }
    explicit IdentityMatrix(int n){
        rows = columns = n;
        data = new double* [rows];
        for (int i = 0; i < n; i++){
            data[i] = new double [rows];
        }
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (i == j){
                    data[i][j] = 1.00;
                } else {
                    data[i][j] = 0.00;
                }
            }
        }
    }
};
class EliminationMatrix: public IdentityMatrix{
public:
    EliminationMatrix(auto &squaredMatrix, int n, int m){
        auto *identityMatrix = new IdentityMatrix(squaredMatrix.getRows());
        identityMatrix->setElement(n, m, -((float)squaredMatrix.getElement(n, m)/(float)squaredMatrix.getElement(m, m)));
        rows = columns = squaredMatrix.getRows();
        data = new double* [rows];
        for (int i = 0; i < rows; i++){
            data[i] = new double [rows];
        }
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                data[i][j] = identityMatrix->getElement(i, j);
            }
        }
    }
};
class PermutationMatrix: public IdentityMatrix{
public:
    PermutationMatrix(auto &squaredMatrix, int n, int m){
        auto *identityMatrix = new IdentityMatrix(squaredMatrix.getRows());
        rows = columns = squaredMatrix.getColumns();
        data = new double*[rows];
        for (int i = 0; i < rows; i++){
            data[i] = new double [rows];
        }
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                if (i == n){
                    if (j == m){
                        data[i][j] = 1;
                    } else {
                        data[i][j] = 0;
                    }
                } else if (i == m){
                    if (j == n){
                        data[i][j] = 1;
                    } else {
                        data[i][j] = 0;
                    }
                } else if (i == j){
                    data [i][j] = 1;
                } else {
                    data[i][j] = 0;
                }
            }
        }
    }
};
class ColumnVector {
private:
    int rows{};
    double *data{};
public:
    ColumnVector(){
        rows = 0;
    }
    explicit ColumnVector(int n){
        rows = n;
        data = new double[rows];
        for (int i = 0; i < rows; i++){
            data[i] = 0.00;
        }
    }
    void setElement(int i, double num){
        data[i] = num;
    }
    friend istream& operator>>(istream&in, ColumnVector &columnVector){
        for(int i = 0; i < columnVector.rows; i++){
            string s;
            in >> s;
            columnVector.data[i] = stod(s);
        }
        return in;
    }
    friend ostream& operator<<(ostream& out, ColumnVector &columnVector) {
        for(int i = 0; i < columnVector.rows; i++){
            out.precision(4);
            out << fixed << columnVector.data[i] << endl;
        }
        return out;
    }
    ColumnVector operator*(squaredMatrix otherMatrix){
        ColumnVector temp(otherMatrix.getRows());
        for (int i = 0; i < this->rows; i++) {
            double t = 0.00;
            for (int j = 0; j < this->rows; j++) {
                t += otherMatrix.getElement(i,j) * this->data[j];
            }
            temp.data[i] = t;
            if (fabs(temp.data[i]) < limit ){
                temp.data[i] = 0.00;
            }
        }
        return temp;
    }

    void normalization(squaredMatrix &otherMatrix){
        for (int i = 0; i < otherMatrix.getRows(); i++) {
            this->data[i] = this->data[i]/otherMatrix.getElement(i,i);
            if (fabs(this->data[i]) < limit ){
                this->data[i] = 0.00;
            }
            otherMatrix.setElement(i, i, 1.00);
        }
    }
};
double calculateE(Matrix &matrix){
    double sum = 0;
    for (int i = 0; i < matrix.getRows(); i++){
        for (int j = 0; j < matrix.getColumns(); j++){
            sum += pow(matrix.getElement(i,j), 2);
        }
    }
    return sqrt(sum);
}
squaredMatrix findInverse(Matrix &matrix){
    int n = matrix.getRows();
    auto *identityMatrix = new IdentityMatrix(n);
    squaredMatrix square = *identityMatrix;
    for (int k = 0; k < n - 1; k++) {
        int f = matrix.findMaxColumn(k);
        if (f != k) {
            PermutationMatrix *permutationMatrix = new PermutationMatrix(matrix, k, f);
            square = *permutationMatrix * square;
            matrix = *permutationMatrix * matrix;
        }
        bool flag = false;
        for (int i = k + 1; i < n; i++){
            if (fabs(matrix.getElement(i, k)) >= limit){
                flag = true;
                break;
            }
        }
        if (flag) {
            for (int i = 0; i < n - k - 1; i++) {
                EliminationMatrix *eliminationMatrix;
                eliminationMatrix = new EliminationMatrix(matrix, k + i + 1, k);
                matrix = *eliminationMatrix * matrix;
                square = *eliminationMatrix * square;
            }
        }
    }
    for (int k = n - 1; k > 0; k--){
        for (int i = k - 1; i >= 0; i--){
            EliminationMatrix *eliminationMatrix;
            eliminationMatrix = new EliminationMatrix(matrix, i, k);
            matrix = *eliminationMatrix * matrix;
            square = *eliminationMatrix * square;
        }
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double delimiter = matrix.getElement(i, i);
            square.setElement(i,j,square.getElement(i,j)/delimiter);
        }
        matrix.setElement(i,i, 1.00);
    }
    return square;
}

/*
 * Input format
The input contains:
the length m of data set
m lines with experimental data ti *space* bi
the degree of the polynomial n

 * Output format
The output contains:
the matrix A itself after the line "A:"
the matrix (AT A) after the line "A_T*A:"
the matrix (AT A)-1 after the line "(A_T*A)_-1:"
the matrix AT b after the line "A_T*b:"
the answer itself after the line "x:"
 */

void least_square_approximation(){
    int m;
    cin >> m;
    vector<int> t;
    vector<int> b;
    for (int i = 0; i < m; i++){
        int ti, bi;
        cin >> ti >> bi;
        t.push_back(ti);
        b.push_back(bi);
    }
    int nn;
    cin >> nn;
    Matrix matrixA(m, nn + 1);
    Matrix columnVector(m, 1);
    for (int i = 0; i < m; i++){
        for (int j = 0; j < nn + 1; j++){
            matrixA.setElement(i, j, pow(t[i], j));
        }
        columnVector.setElement(i, 0, b[i]);
    }
    cout << "A:\n" << matrixA;
    Matrix At = matrixA.transpose(matrixA);
    Matrix mulTransA = At * matrixA;
    cout << "A_T*A:\n" << mulTransA;
    squaredMatrix square = findInverse(mulTransA);
    cout << "(A_T*A)^-1:\n" << square;
    Matrix mulTransCol = At * columnVector;
    cout << "A_T*b:\n" << mulTransCol ;
    Matrix x = square * mulTransCol;
    cout << "x~:\n" << x;
}

/*
 * Input format
The input contains:
A square matrix A (in element-wise manner with the dimension firstly) as in the previous exercises.
A vector of free coefficients b (in element-wise manner with the dimension firstly).
The approximation accuracy.

 * Output format
The output contains:
The string "The method is not applicable!"
or
Matrix alpha, entitled "alpha:"
Vector beta, entitled "beta:"
Set of vectors xi of the approximation steps, each entitled "x(i):"
Current accuracy e for each step, entitled "e:" (skipped for the last step).
 */

void jacobi_method(){
    int n;
    cin >> n;
    squaredMatrix A (n);
    cin >> A;
    cin >> n;
    Matrix B(n, 1);
    cin >> B;
    double e;
    cin >> e;
    squaredMatrix D(n);
    for (int i = 0; i < n; i++){
        D.setElement(i, i, A.getElement(i, i));
    }
    bool flag = true;
    for (int i = 0; i < n; i++){
        double sum = 0;
        for (int j = 0; j < n; j++){
            if (i != j){
                sum += fabs(A.getElement(i, j));
            }
        }
        if (fabs(A.getElement(i,i)) < sum){
            flag = false;
            break;
        }
    }
    if (!flag){
        cout << "The method is not applicable!";
    } else {
        squaredMatrix D_I(n);
        for (int i = 0; i < n; i++) {
            D_I.setElement(i, i, 1 / D.getElement(i, i));
        }
        squaredMatrix alpha(n);
        IdentityMatrix I(n);
        squaredMatrix multiply = D_I * A;
        alpha = I - multiply;
        Matrix beta(n, 1);
        cout << "alpha:\n" << alpha;
        beta = D_I * B;
        cout << "beta:\n" << beta;
        Matrix x0 = beta;
        cout << "x(0):\n" << x0;
        squaredMatrix difference = (A - D) * x0;
        double ee = 1;
        int iteration = 0;
        while (ee > e) {
            iteration++;
            Matrix xNext = D_I * (B - difference);
            difference = (A - D) * xNext;
            Matrix deltaX = xNext - x0;
            ee = calculateE(deltaX);
            cout << "e: " << ee << endl;
            x0 = xNext;
            cout << "x(" << iteration <<"):\n" << xNext;
        }

    }
};

/*
 * Input format
The input contains:

A square matrix A (in element-wise manner with the dimension firstly) as in the previous exercises.
A vector of free coefficients b (in element-wise manner with the dimension firstly).
The approximation accuracy .
 * Output format
The output contains:
The string "The method is not applicable!"
or
Matrix alpha, entitled "alpha:"
Vector beta, entitled "beta:"
Matrix B, entitled "B:"
Matrix C, entitled "C:"
Matrix I−B, entitled "I-B:"
Matrix (I−B)^-1, entitled "(I-B)_-1:"
 */

void seidel_method(){
    int n;
    cin >> n;
    squaredMatrix A (n);
    cin >> A;
    cin >> n;
    Matrix b(n, 1);
    cin >> b;
    double e;
    cin >> e;
    squaredMatrix D(n);
    for (int i = 0; i < n; i++){
        D.setElement(i, i, A.getElement(i, i));
    }
    bool flag = true;
    for (int i = 0; i < n; i++){
        double sum = 0;
        for (int j = 0; j < n; j++){
            if (i != j){
                sum += fabs(A.getElement(i, j));
            }
        }
        if (fabs(A.getElement(i,i)) < sum){
            flag = false;
            break;
        }
    }
    if (!flag){
        cout << "The method is not applicable!";
    } else {
        squaredMatrix D_I(n);
        for (int i = 0; i < n; i++) {
            D_I.setElement(i, i, 1 / D.getElement(i, i));
        }
        squaredMatrix alpha(n);
        IdentityMatrix I(n);
        squaredMatrix multiply = D_I * A;
        alpha = I - multiply;
        squaredMatrix C(n);
        squaredMatrix B(n);
        squaredMatrix L(n), U(n);
        for (int j = 0; j < n ; j++){
            for (int i = 0; i < n; i++){
                if (i > j){
                    C.setElement(j , i, alpha.getElement(j, i));
                    U.setElement(j , i, A.getElement(j, i));
                }
            }
        }
        for (int j = 0; j < n ; j++){
            for (int i = 0; i < n; i++){
                if (i < j){
                    B.setElement(j , i, alpha.getElement(j, i));
                }
                if (i<=j){
                    L.setElement(j , i, A.getElement(j, i));
                }
            }
        }
        Matrix beta(n, 1);
        beta = D_I * b;
        cout << "beta:\n" << beta;
        Matrix x0 = beta;
        cout << "alpha:\n" << alpha;
        cout << "B:\n" << B;
        cout << "C:\n" << C;
        squaredMatrix I_B = I - B;
        cout << "I-B:\n" << I_B;
        squaredMatrix I_B_1 = findInverse(I_B);
        cout << "(I-B)_-1:\n" << I_B_1;
        cout << "x(0):\n" << x0;
        Matrix difference = b - (U * x0);
        squaredMatrix LI = findInverse(L);
        double ee = 1;
        int iteration = 0;
        while (ee > e) {
            iteration++;
            Matrix xNext = LI * difference;
            difference = b - (U * xNext);
            Matrix deltaX = xNext - x0;
            ee = calculateE(deltaX);
            cout << "e: " << ee << endl;
            x0 = xNext;
            cout << "x(" << iteration <<"):\n" << xNext;
        }

    }
}

int main() {
    return 0;
}