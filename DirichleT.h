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
	TMatrix<double> V; //������� ������
	TVector<double> xBorder, yBorder; //������� �� ��� x � y ��������������
	double h, k; //��� �� ��� x � y ��������������
	int TaskNumber, n, m; //����� ������, ����������� �����
	double w; //�������� ������
public:
	DirichleT(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber); //����������� �������������

	double F_Fuction(double x, double y); //������� � ������ ����� (��� ������, ������ f(x,y))
	double ExactSolution(double x, double y);//������ �������
	double XInicialConditions(double x, int Num); //��������� ������� � ���� ������� ����������� ��� X
	double YInicialConditions(double y, int Num); //��������� ������� � ���� ������� ����������� ��� Y
	void Inicialisation(); //���������� ��������� ������� � ���� �������

	double UpRelaxMethod(); //����� ������� ����������

	TVector<double> SetOptW(); //���������� ����������� ���
	void SetW(double W); //���������� ���������������� ���

	void GetRes(double **mas); //�������� ������� � ������������ ������. �����, ����� ���������� ����� � ������

	TVector<double> MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations); //�������� �������
	TVector<double> MethodError(double eps, int MaxIterations); //����������� �������

	double NevyazkaInf(); //������� �� ����� �������������
	double NevyazkaEvkl(); //��������� ����� �������

	void ChangeGrid(int N, int M); //�������� �����
	void ChangeBorders(TVector<double> XBorder, TVector<double> YBorder); //�������� �������
	void ChangeTask(int TASKNumber); //�������� ������
	void Reset();//��������� �����������

};

#endif