//
// Created by Daniel on 07.05.2019.
//
#include <iostream>
#include "SquareMatrix.h"
#include "cdecomp.h"
#include <cmath>

SquareMatrix::~SquareMatrix() {
    delete[] matrix;
}

/*Создает матрицу, проверяя ее порядок на правильность при этом заполняет главную диогональ нулями*/
SquareMatrix::SquareMatrix(int order)
        : order(order), matrix(new double[order * order]), cond(0), pivot(new int[order]) {
    if (order < 1) {
        throw std::invalid_argument("SquareMatix::Constructor: wrong order of unit matrix");
    }
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            this->matrix[i * order + j] = (i == j);
        }
    }
}

/*Создет матрицу, используя данный порядок и считывая матрицу из массива. Массив задает матрицу построчно*/
SquareMatrix::SquareMatrix(int order, const double *matrix)
        : order(order), matrix(new double[order * order]), cond(0), pivot(new int[order]) {
    if (order < 1) {
        throw std::invalid_argument("SquareMatix::Constructor: wrong order of matrix");
    }
    for (int i = 0; i < order * order; ++i) {
        this->matrix[i] = matrix[i];
    }
}

SquareMatrix::SquareMatrix(int order, std::ifstream& in)
        : order(order), matrix(new double[order * order]), cond(0), pivot(new int[order]) {
    if (order < 1) {
        throw std::invalid_argument("SquareMatix::Constructor: wrong order of matrix");
    }
    for (int i = 0; i < order * order; i++) {
        in >> matrix[i];
    }
}

/*Конструктор копирования*/
SquareMatrix::SquareMatrix(SquareMatrix const &squareMatrix)
        : order(squareMatrix.order), matrix(new double[order * order]), cond(squareMatrix.cond), pivot(new int[order]) {
    for (int i = 0; i < order; ++i) {
        pivot[i] = squareMatrix.pivot[i];
    }
    for (int i = 0; i < order * order; ++i) {
        this->matrix[i] = squareMatrix.matrix[i];
    }
}

/*Перегрузка оператора присваивания*/
SquareMatrix &SquareMatrix::operator=(SquareMatrix const &squareMatrix) {
    if (this != &squareMatrix) {
        order = squareMatrix.order;
        cond = squareMatrix.cond;
        delete[] pivot;
        pivot = new int[order];
        for (int i = 0; i < order; ++i) {
            pivot[i] = squareMatrix.pivot[i];
        }
        delete[] matrix;
        matrix = new double[order * order];
        for (int i = 0; i < order * order; ++i) {
            matrix[i] = squareMatrix.matrix[i];
        }
    }
    return *this;
}

/*Константные геттеры для матрицы и порядка*/
double *SquareMatrix::getMatrix() const {
    return matrix;
}

int SquareMatrix::getOrder() const {
    return order;
}

double SquareMatrix::getCond() const {
    return cond;
}

int *SquareMatrix::getPivot() const {
    return pivot;
}

SquareMatrix &SquareMatrix::operator+=(SquareMatrix const &squareMatrix) {
    if (order != squareMatrix.order) {
        throw std::invalid_argument("SquareMatix::operator+= :wrong orders in sum");
    }
    for (int i = 0; i < order * order; ++i) {
        matrix[i] += squareMatrix.matrix[i];
    }
    return *this;
}

SquareMatrix &SquareMatrix::operator-=(SquareMatrix const &squareMatrix) {
    if (order != squareMatrix.order) {
        throw std::invalid_argument("SquareMatix::operator-=: wrong matrix orders");
    }
    for (int i = 0; i < order * order; ++i) {
        matrix[i] -= squareMatrix.matrix[i];
    }
    return *this;
}

SquareMatrix &SquareMatrix::operator*=(SquareMatrix const &squareMatrix) {
    SquareMatrix newSquareMatrix = *this;
    newSquareMatrix = newSquareMatrix * squareMatrix;
    *this = newSquareMatrix;
    return *this;
}

void SquareMatrix::printMatrix() {
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            if (matrix[i * order + j] < 0) {
                std::cout << matrix[i * order + j] << " ";
            } else {
                std::cout << " " << matrix[i * order + j] << " ";
            }
        }
        std::cout << std::endl;
    }
}

void SquareMatrix::convertToLU() {
    int flag = 0;
    decomp(order, order, matrix, &cond, pivot, &flag);

    if (flag == 1) {
        throw std::invalid_argument("SquareMatix::convertToLU: could not allocate memory for workspace.");
    }
    if (flag == 2) {
        throw std::invalid_argument(
                "SquareMatix::convertToLU: illegal user input n < 1, a == NULL, pivot == NULL, n > ndim.");
    }
    if (flag == 3) {
        throw std::invalid_argument("SquareMatix::convertToLU: matrix is singular");
    }
}

double SquareMatrix::operator()(int i, int j) {
    if ((i >= order) || (i < 0) || (j >= order) || (j < 0)) {
        throw std::out_of_range("SquareMatix::operator(): i or j out of range");
    }
    return matrix[i * order + j];
}

SquareMatrix &SquareMatrix::operator!() {
    SquareMatrix temp = SquareMatrix(*this);
    temp.convertToLU();
    cond = temp.cond;
    pivot = temp.pivot;
    auto *result = new SquareMatrix(order);
    auto *temp2 = new double[order];
    for (int i = 0; i < order; ++i) {
        result->getColumn(i, temp2);
        solve(order, order, temp.getMatrix(), temp2, temp.getPivot());
        result->setColumn(i, temp2);
    }
    delete[] temp2;
    return *result;
}

SquareMatrix &SquareMatrix::operator~() {
    auto *result = new SquareMatrix(order);
    auto *temp = new double[order];
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            temp[j] = matrix[i * order + j];
        }
        result->setColumn(i, temp);
    }
    delete[] temp;
    return *result;
}

void SquareMatrix::getColumn(int i, double *res) {
    if ((i >= order) || (i < 0)) {
        throw std::out_of_range("SquareMatix::operator(): i out of range");
    }
    for (int j = 0; j < order; ++j) {
        res[j] = matrix[(j * order) + i];
    }
}

void SquareMatrix::setColumn(int i, const double *source) {
    if ((i >= order) || (i < 0)) {
        throw std::out_of_range("SquareMatix::operator(): i out of range");
    }
    for (int j = 0; j < order; ++j) {
        matrix[(j * order) + i] = source[j];
    }
}

double SquareMatrix::getRate() {
    double sum = 0;
    for (int i = 0; i < order * order; ++i) {
        sum += (matrix[i] * matrix[i]);
    }
    return sqrt(sum);
}

void SquareMatrix::Givens_method() {
    for (int j = 0; j < order - 2; j++) {
        for (int i = j + 2; i < order; i++) {
            SquareMatrix T = getTij(i, j + 1);
            *this = ~T * (*this * T);
        }
    }
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            if (i > j + 1) {
                matrix[i * order + j] = abs(round(matrix[i * order + j] * 100000000000000) / 100000000000000);
            }
        }
    }
}

SquareMatrix SquareMatrix::getTij(int i, int j) {
    double C = matrix[j * order + (j - 1)] / (sqrt(matrix[j * order + (j - 1)] * matrix[j * order + (j - 1)] +
                                                   matrix[i * order + (j - 1)] * matrix[i * order + (j - 1)]));
    double S = -(matrix[i * order + (j - 1)]) / (sqrt(matrix[j * order + (j - 1)] * matrix[j * order + (j - 1)] +
                                                      matrix[i * order + (j - 1)] * matrix[i * order + (j - 1)]));
    SquareMatrix temp(order);
    temp.setElem(i, i, C);
    temp.setElem(j, j, C);
    temp.setElem(j, i, S);
    temp.setElem(i, j, -S);
    return temp;
}

void SquareMatrix::setElem(int i, int j, double source) {
    if (i >= order || j >= order || i < 0 || j < 0) {
        throw std::invalid_argument("SquareMatrix::setElem: wrong i or j");
    }
    matrix[i * order + j] = source;
}

SquareMatrix operator+(SquareMatrix const &squareMatrix1, SquareMatrix const &squareMatrix2) {
    SquareMatrix newSquareMatrix = squareMatrix1;
    newSquareMatrix += squareMatrix2;
    return newSquareMatrix;
}

SquareMatrix operator-(SquareMatrix const &squareMatrix1, SquareMatrix const &squareMatrix2) {
    SquareMatrix newSquareMatrix = squareMatrix1;
    newSquareMatrix -= squareMatrix2;
    return newSquareMatrix;
}

SquareMatrix operator*(SquareMatrix const &squareMatrix1, SquareMatrix const &squareMatrix2) {
    if (squareMatrix1.getOrder() != squareMatrix2.getOrder()) {
        throw std::invalid_argument("SquareMatix::operator*: wrong matrix orders");
    }
    /*Получаем указатели на массивы, дабы избежать множественных вызовов*/
    int order = squareMatrix1.getOrder();
    auto *arr1 = squareMatrix1.getMatrix();
    auto *arr2 = squareMatrix2.getMatrix();
    auto *res = new double[order * order];
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            double sum = 0;
            for (int k = 0; k < order; ++k) {
                sum += arr1[(i * order) + k] * arr2[(k * order) + j];
            }
            res[(i * order) + j] = sum;
        }
    }
    SquareMatrix newSquareMatrix(order, res);
    delete[] res;
    return newSquareMatrix;
}