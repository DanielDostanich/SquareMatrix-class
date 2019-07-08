//
// Created by Daniel on 07.05.2019.
//

#ifndef __SQUARE_MATRIX_H
#define __SQUARE_MATRIX_H

#include <fstream>

class SquareMatrix {
private:
    /*Матрица хранится как одномерный массив. Размер строки задается порядком.*/
    int order;
    double * matrix;
    /*Как поля хранятся данные необходимые для работы функции solve (pivot), а также число обусловленности матрицы.
     * До вызова convertToLU либо operator!() хранят в себе неопределенные данные*/
    double cond;
    int * pivot;
public:
    /*Деструктор*/
    ~SquareMatrix();
    /*Дефолтный конструктор*/
    SquareMatrix() = default;
    /*Создание единичной матрицы*/
    explicit SquareMatrix(int order);
    /*Создание матрицы из массива*/
    SquareMatrix(int order, const double * matrix);
    /*Конструктор копирования*/
    SquareMatrix (SquareMatrix const & squareMatrix);
    SquareMatrix (int order, std::ifstream& in);
    /*Оператор присваивания*/
    SquareMatrix &operator=(SquareMatrix const& squareMatrix);
    /*Перегруженные операторы для сложения/вычитания/умножения*/
    SquareMatrix &operator+=(SquareMatrix const& squareMatrix);
    SquareMatrix &operator-=(SquareMatrix const& squareMatrix);
    SquareMatrix &operator*=(SquareMatrix const& squareMatrix);
    /*Получение из матрицы ее обратной*/
    SquareMatrix &operator!();
    /*Получение из матрицы транспонированной*/
    SquareMatrix &operator~();
    double operator()(int i, int j);
    /*Получение массива элементов из матрицы*/
    double * getMatrix() const;
    /*Получение порядка матрицы*/
    int getOrder() const;
    /*Получение cond и pivot. До использования convertToLU или operator!() хранят в себе неопределенные данные*/
    double getCond() const;
    int * getPivot() const;
    SquareMatrix getTij(int i, int j);
    /*Помещает значения столбца i в массив res*/
    void getColumn(int i, double * res);
    /*Помещает значения из массива source в столбец i*/
    void setColumn(int i, double const * source);
    /*Распечатывание матрицы*/
    void printMatrix();
    /*Приводит матрицу к форме LU разложения*/
    void convertToLU();
    double getRate();
    void setElem(int i, int j, double source);
    void Givens_method();
};
/*Перегруженные операторы для сложения/вычитания/умножения*/
SquareMatrix operator+(SquareMatrix const& squareMatrix1, SquareMatrix const& squareMatrix2);
SquareMatrix operator-(SquareMatrix const& squareMatrix1, SquareMatrix const& squareMatrix2);
SquareMatrix operator*(SquareMatrix const& squareMatrix1, SquareMatrix const& squareMatrix2);

#endif //__SQUARE_MATRIX_H