#include <iostream>
#include <utility>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

class Matrix{
protected:
    int rows;
    int columns;
    vector <vector <double>> elements;
public:
    Matrix() {
        rows = 0;
        columns = 0;
    }
    Matrix (int newRows, int newColumns, vector <vector <double>> newElements) {
        rows = newRows;
        columns = newColumns;
        elements = std::move(newElements);
    }
    double get(int row, int column) {
        return elements[row][column];
    }
    int getRows() const {
        return rows;
    }
    int getColumns() const {
        return columns;
    }
    vector <vector <double>> getElements(){
        return elements;
    }
    Matrix& operator= (Matrix move) {
        if (this == &move) {
            return *this;
        }
        rows = move.getRows();
        columns = move.getColumns();
        elements = move.getElements();
        return *this;
    }
    Matrix operator+ (Matrix sum) {
        if (rows != sum.getRows() || columns!= sum.getColumns()) {
            throw std::invalid_argument("Error: the dimensional problem occurred\n");
        }
        vector <vector <double>> newElements(rows);
        for (int i = 0; i < rows; i++) {
            newElements[i].resize(columns);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                newElements[i][j] = elements[i][j] + sum.getElements()[i][j];

            }
        }
        Matrix result(rows, columns, newElements);
        return result;
    }
    Matrix operator- (Matrix sub) {
        if (rows != sub.getRows() || columns != sub.getColumns()) {
            throw std::invalid_argument("Error: the dimensional problem occurred\n");
        }
        vector <vector <double>> newElements(rows);
        for (int i = 0; i < rows; i++) {
            newElements[i].resize(columns);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                newElements[i][j] = elements[i][j] - sub.getElements()[i][j];
            }
        }
        Matrix result(rows, columns, newElements);
        return result;
    }

    virtual Matrix operator* (Matrix mul) {
        if (columns != mul.getRows()) {
            throw std::invalid_argument("Error: the dimensional problem occurred\n");
        }
        vector <vector <double>> newElements(rows);
        for (int i = 0; i < rows; i++) {
            newElements[i].resize(mul.getColumns());
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < mul.getColumns(); j++) {
                double element = 0;
                for (int k = 0; k < columns; k++) {
                    element += elements[i][k] * mul.getElements()[k][j];
                }
                newElements[i][j] = element;
            }
        }
        Matrix result(rows, mul.getColumns(), newElements);
        return result;
    }
    Matrix transpose() {
        vector <vector <double>> newElements(columns);
        for (int i = 0; i < columns; i++) {
            newElements[i].resize(rows);
        }
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                newElements[i][j] = elements[j][i];
            }
        }
        Matrix result(columns, rows, newElements);
        return result;
    }
    friend ostream& operator<<(ostream& os, Matrix matrix) {
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                os << matrix.elements[i][j];
                if (j != matrix.columns - 1) {
                    os << " ";
                }
            }
            os << endl;
        }
        return os;
    }

    friend istream& operator>>(istream& is, Matrix &matrix) {
        is >> matrix.rows >> matrix.columns;
        matrix.elements.resize(matrix.rows);
        for (int i = 0; i < matrix.rows; i++) {
            matrix.elements[i].resize(matrix.columns);
        }
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                is >> matrix.elements[i][j];
            }
        }
        return is;
    }
};

class SquareMatrix: public Matrix {
public:
    SquareMatrix() : Matrix() {}

    SquareMatrix(int size, vector <vector <double>> newElements) : Matrix(size, size,std::move(newElements)) {}

    SquareMatrix& operator= (Matrix move) {
        if (this == &move) {
            return *this;
        }
        rows = move.getRows();
        columns = move.getColumns();
        elements = move.getElements();
        return *this;
    }

    friend istream& operator>>(istream& is, SquareMatrix &matrix) {
        is >> matrix.rows;
        matrix.columns = matrix.rows;
        matrix.elements.resize(matrix.rows);
        for (int i = 0; i < matrix.rows; i++) {
            matrix.elements[i].resize(matrix.columns);
        }
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.columns; j++) {
                is >> matrix.elements[i][j];
            }
        }
        return is;
    }
};

class IdentityMatrix: public SquareMatrix {
public:
    IdentityMatrix() : SquareMatrix() {}

    explicit IdentityMatrix(int size) : SquareMatrix() {
        rows = size;
        columns = size;
        elements = vector<vector<double>>(rows);
        for (int i = 0; i < rows; i++) {
            elements[i].resize(columns);
            for (int j = 0; j < rows; j++) {
                elements[i][j] = i == j ? 1 : 0;
            }
        }
    }
};

class EliminationMatrix : public SquareMatrix {
private:
    int nullifyRow;
    int nullifyColumn;
public:
    EliminationMatrix() : SquareMatrix() {
        nullifyRow = 0;
        nullifyColumn = 0;
    }

    EliminationMatrix(int rowNumber, int columnNumber, Matrix &toNullify) : SquareMatrix() {
        rows = toNullify.getRows();
        columns = toNullify.getRows();
        nullifyRow = rowNumber - 1;
        nullifyColumn = columnNumber - 1;
        elements = vector<vector<double>>(rows);
        for (int i = 0; i < rows; i++) {
            elements[i].resize(columns);
            for (int j = 0; j < rows; j++) {
                elements[i][j] = i == j ? 1 : 0;
            }
        }
        elements[nullifyRow][nullifyColumn] = - toNullify.get(nullifyRow, nullifyColumn)
                                              / toNullify.get(nullifyColumn, nullifyColumn);
    }
};

class PermutationMatrix : public SquareMatrix {
private:
    int permuteRow1;
    int permuteRow2;
public:
    PermutationMatrix() : SquareMatrix() {
        permuteRow1 = 0;
        permuteRow2 = 0;
    }
    PermutationMatrix(int row1, int row2, Matrix &toPermute) : SquareMatrix() {
        rows = toPermute.getRows();
        columns = toPermute.getRows();
        permuteRow1 = row1 - 1;
        permuteRow2 = row2 - 1;
        elements = vector<vector<double>>(rows);
        for (int i = 0; i < rows; i++) {
            elements[i].resize(columns);
            for (int j = 0; j < rows; j++) {
                elements[i][j] = i == j ? 1 : 0;
            }
        }
        elements[permuteRow1][permuteRow2] = 1;
        elements[permuteRow2][permuteRow1] = 1;
        elements[permuteRow1][permuteRow1] = 0;
        elements[permuteRow2][permuteRow2] = 0;
    }
};

class NormalizationMatrix : public SquareMatrix {
public:
    NormalizationMatrix() : SquareMatrix() {}
    explicit NormalizationMatrix(Matrix &toNormalize) : SquareMatrix() {
        rows = toNormalize.getRows();
        columns = rows;
        elements = vector<vector<double>>(rows);
        for (int i = 0; i < rows; i++) {
            elements[i].resize(columns);
            for (int j = 0; j < rows; j++) {
                elements[i][j] = i == j ? (1.00 / toNormalize.get(i, i)) : 0;
            }
        }
    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector() : Matrix() {}
    ColumnVector(int size, vector <vector <double>> newElements) : Matrix(size, 1,std::move(newElements)) {}
    double norm() {
        Matrix N = *this * *this;
        return N.get(0, 0);
    }
    friend istream& operator>>(istream& is, ColumnVector &vector) {
        is >> vector.rows;
        vector.columns = 1;
        vector.elements.resize(vector.rows);
        for (int i = 0; i < vector.rows; i++) {
            vector.elements[i].resize(vector.columns);
        }
        for (int i = 0; i < vector.rows; i++) {
            for (int j = 0; j < vector.columns; j++) {
                is >> vector.elements[i][j];
            }
        }
        return is;
    }
    ColumnVector& operator= (Matrix move) {
        if (this == &move) {
            return *this;
        }
        rows = move.getRows();
        columns = move.getColumns();
        elements = move.getElements();
        return *this;
    }
};

int main() {
    cout << fixed;
    cout << setprecision(4);
    int pointsNumber;
    cin >> pointsNumber;
    vector<double> xValues{};
    vector<double> yValues{};
    for (int i = 0; i < pointsNumber; i++) {
        int x, y;
        cin >> x >> y;
        xValues.push_back(x);
        yValues.push_back(y);
    }
    int size;
    cin >> size;
    vector<vector<double>> elements;
    vector<vector<double>> vectorElements;
    vectorElements.resize(pointsNumber);
    elements.resize(pointsNumber);
    for (int i = 0; i < pointsNumber; i++) {
        elements[i].resize(size + 1);
        vectorElements[i].resize(1);
        vectorElements[i][0] = yValues[i];
        for (int j = 0; j < size + 1; j++) {
            elements[i][j] = pow(xValues[i], j);
        }
    }
    Matrix A{pointsNumber, size + 1, elements};
    ColumnVector b{pointsNumber, vectorElements};
    cout << "A:" << endl;
    cout << A;
    Matrix ATransposed = A.transpose();
    Matrix ATransposedA = ATransposed * A;
    cout << "A_T*A:" << endl;
    cout << ATransposedA;
    Matrix ATransposedAInverse = IdentityMatrix(ATransposedA.getRows());
    for (int i = 0; i < ATransposedA.getRows(); i++) {
        double maxPivot = abs(ATransposedA.get(i, i));
        int maxPivotRow = i;
        for (int j = i; j < ATransposedA.getColumns(); j++) {
            if (maxPivot < abs(ATransposedA.get(j, i))) {
                maxPivot = abs(ATransposedA.get(j, i));
                maxPivotRow = j;
            }
        }
        if (maxPivotRow != i) {
            PermutationMatrix P(i + 1, maxPivotRow + 1, ATransposedA);
            ATransposedA = P * ATransposedA;
            ATransposedAInverse = P * ATransposedAInverse;
        }
        for (int j = i + 1; j < ATransposedA.getColumns(); j++) {
            if (ATransposedA.get(j, i) != 0) {
                EliminationMatrix E(j + 1, i + 1, ATransposedA);
                ATransposedA = E * ATransposedA;
                ATransposedAInverse = E * ATransposedAInverse;
            }
        }
    }
    for (int i = ATransposedA.getColumns() - 1; i > -1; i--) {
        for (int j = i - 1; j > -1; j--) {
            if (ATransposedA.get(j, i) != 0) {
                EliminationMatrix E(j + 1, i + 1, ATransposedA);
                ATransposedA = E * ATransposedA;
                ATransposedAInverse = E * ATransposedAInverse;
            }
        }
    }
    NormalizationMatrix N(ATransposedA);
    ATransposedA = N * ATransposedA;
    ATransposedAInverse = (N * ATransposedAInverse);
    cout << "(A_T*A)^-1:" << endl;
    cout << ATransposedAInverse;
    Matrix ATransposedB = ATransposed * b;
    cout << "A_T*b:" << endl;
    cout << ATransposedB;
    Matrix x = ATransposedAInverse * ATransposedB;
    cout << "x~:" << endl;
    cout << x;
    return 0;
}

