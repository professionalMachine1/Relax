#pragma once

#ifndef Dirichle_T
#define DirichleT_H

#include <iostream>
#include <math.h>
#include "TMatrix.h"

using namespace std;

class DirichleT
{
private:
	TMatrix<double> V; //Решение задачи
	TVector<double> xBorder, yBorder; //Границы по оси x и y соответственно
	double h, k; //Шаг по оси x и y соответственно
	int TaskNumber, n, m; //Номер задачи, размерность сетки
	double w; //Параметр метода
public:
	DirichleT(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber); //Конструктор инициализации

	double F_Fuction(double x, double y); //Функция в правой части (без минуса, только f(x,y))
	double ExactSolution(double x, double y);//Точные решения
	double XInicialConditions(double x, int Num); //Начальные условия в виде функций параллельно оси X
	double YInicialConditions(double y, int Num); //Начальные условия в виде функций параллельно оси Y
	void Inicialisation(); //Заполнение начальных условий в наше решение

	double UpRelaxMethod(); //Метод верхней релаксации

	TVector<double> SetOptW(); //Установить оптимальный шаг
	void SetW(double W); //Установить пользовательский шаг

	void GetRes(double **mas); //Записать решение в динамический массив. Нужно, чтобы съэконоить время и память

	TVector<double> MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations); //Точность решения
	TVector<double> MethodError(double eps, int MaxIterations); //Погрешность решения

	double NevyazkaInf(); //Невязка по норме бесконечность
	double NevyazkaEvkl(); //Евклидова норма невязки

	void ChangeGrid(int N, int M); //Изменить сетку
	void ChangeBorders(TVector<double> XBorder, TVector<double> YBorder); //Изменить границу
	void ChangeTask(int TASKNumber); //Изменить задачу
	void Reset();//Обнуление результатов

};

#endif