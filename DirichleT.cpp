#include "DirichleT.h"


double pi = 2 * asin(1);


DirichleT::DirichleT(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber)
{
	TaskNumber = TASKNumber;
	n = N;
	m = M;
	xBorder = XBorder;
	yBorder = YBorder;
	h = (xBorder[xBorder.Size() - 1] - xBorder[0]) / (1.0*n);
	k = (yBorder[yBorder.Size() - 1] - yBorder[0]) / (1.0*m);
	V = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
		V[j] = TVector<double>(n + 1);
	Inicialisation();
	w = 1;
}


double DirichleT::F_Fuction(double x, double y)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		res = pow(pi, 2)*sin(pi*x*y)*(pow(x, 2) + pow(y, 2));
		break;
	case 2:
		res = 4;
		break;
	case 3:
		res = -exp(x*y)*(pow(x, 2) + pow(y, 2));
		break;
	case 4:
		res = -exp(-x * pow(y, 2));
		break;
	default:
		res = 0;
	}
	return res;
}


double DirichleT::ExactSolution(double x, double y)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		res = sin(pi*x*y);
		break;
	case 2:
		res = 1 - pow(x - 1, 2) - pow(y - 0.5, 2);
		break;
	case 3:
		res = exp(x*y);
		break;
	default:
		res = 0;
	}
	return res;
}


double DirichleT::XInicialConditions(double x, int Num)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		switch (Num)
		{
		case 1:
			res = sin(yBorder[0] * pi*x);
			break;
		case 2:
			res = sin(yBorder[yBorder.Size() - 1] * pi*x);
			break;
		}
		break;
	case 2:
		switch (Num)
		{
		case 1:
			res = 0.75 - pow(x - 1, 2);
			break;
		case 2:
			res = 0.75 - pow(x - 1, 2);
			break;
		}
		break;
	case 3:
		switch (Num)
		{
		case 1:
			res = 1;
			break;
		case 2:
			res = exp(x);
			break;
		}
		break;
	case 4:
		switch (Num)
		{
		case 1:
			res = (x - 1)*(x - 2);
			break;
		case 2:
			res = x * (x - 1)*(x - 2);
			break;
		}
		break;
	default:
		res = 0;
	}
	return res;
}


double DirichleT::YInicialConditions(double y, int Num)
{
	double res;
	switch (TaskNumber)
	{
	case 1:
		switch (Num)
		{
		case 1:
			res = sin(xBorder[0] * pi*y);
			break;
		case 2:
			res = sin(xBorder[xBorder.Size() - 1] * pi*y);
			break;
		}
		break;
	case 2:
		switch (Num)
		{
		case 1:
			res = -pow(y - 0.5, 2);
			break;
		case 2:
			res = -pow(y - 0.5, 2);
			break;
		}
		break;
	case 3:
		switch (Num)
		{
		case 1:
			res = exp(-y);
			break;
		case 2:
			res = 1;
			break;
		}
		break;
	case 4:
		switch (Num)
		{
		case 1:
			res = (y - 2)*(y - 3);
			break;
		case 2:
			res = y * (y - 2)*(y - 3);
			break;
		}
		break;
	default:
		res = 0;
	}
	return res;
}


void DirichleT::Inicialisation()
{
	double x, y = yBorder[0];
	for (int j = 0; j <= m; j++)
	{
		V[j][0] = YInicialConditions(y, 1);
		V[j][n] = YInicialConditions(y, 2);
		y += k;
	}
	x = xBorder[0] + h;
	for (int i = 1; i < n; i++)
	{
		V[0][i] = XInicialConditions(x, 1);
		V[m][i] = XInicialConditions(x, 2);
		x += h;
	}
}


TVector<double> DirichleT::SetOptW()
{
	double minSC;
	double temp1, temp2;
	temp1 = pow(sin(pi / (2.0*n)), 2);
	temp2 = pow(sin(pi / (2.0*m)), 2);

	minSC = 4 * temp1 / pow(h, 2) + 4 * temp2 / pow(k, 2); //Минимальное собственное число матрицы A

	temp2 = 2 * pow(k, 2)*temp1 + 2 * pow(h, 2)*temp2;
	temp1 = temp2 / (pow(h, 2) + 2 * pow(k, 2)*temp1);
	w = 2 / (1 + sqrt(temp1*(2 - temp1))); //Оптимальный параметр

	TVector<double> result(2);
	result[0] = w; //Оптмальный параметр
	result[1] = minSC; //Минимальное собственное число матрицы A
	return result;
}


void DirichleT::SetW(double W)
{
	w = W;
}


double DirichleT::UpRelaxMethod()
{
	double A, x, y = yBorder[0] + k, hE, kE;
	double errorM = 0; //Точность метода
	double temp, prev;
	hE = 1 / pow(h, 2);
	kE = 1 / pow(k, 2);
	A = -2 * (hE + kE);
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			prev = V[j][i];
			//
			temp = A * prev + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]);
			V[j][i] = prev - w * (temp + F_Fuction(x, y)) / A;
			//
			prev = fabs(V[j][i] - prev);
			if (prev > errorM)
				errorM = prev;
			x += h;
		}
		y += k;
	}
	return errorM;
}


void DirichleT::GetRes(double **mas)
{
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
			mas[j][i] = V[j][i];
}


TVector<double> DirichleT::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations) //Точность решения
{
	DirichleT Solution(2 * n, 2 * m, xBorder, yBorder, TaskNumber);
	SetOptW();
	Solution.SetOptW();

	TVector<double> accurancy(2); //Точность метода
	//Чтобы зайти в цикл
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];
	TVector<int> IterationsCount(2);
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		accurancy[0] = UpRelaxMethod();
		IterationsCount[0]++;
	}
	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		accurancy[1] = Solution.UpRelaxMethod();
		IterationsCount[1]++;
	}

	double error = 0; //Точность решения
	double sup; //Вспомогательная переменная
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - Solution.V[2 * j][2 * i]);
			if (sup > error)
				error = sup;
		}

	sup = NevyazkaEvkl();
	TVector<double> result(6);
	result[0] = accurancy[0]; //Точность метода на сетке (n+1,m+1)
	result[1] = IterationsCount[0]; //Количество итераций на сетке (n+1,m+1)
	result[2] = accurancy[1]; //Точность метода на сетке (2n+1,2m+1)
	result[3] = IterationsCount[1]; //Количество итераций на сетке (2n+1,2m+1)
	result[4] = error; //Точность решения
	result[5] = sup; //Евклидова норма невязки
	
	return result;
}


TVector<double> DirichleT::MethodError(double eps, int MaxIterations) //Погрешность решения
{
	//Считаем точное решение
	TMatrix<double> u(m + 1);
	double x, y = yBorder[0];
	for (int j = 0; j <= m; j++)
	{
		u[j] = TVector<double>(n + 1);
		x = xBorder[0];
		for (int i = 0; i <= n; i++)
		{
			u[j][i] = ExactSolution(x, y);
			x += h;
		}
		y += k;
	}

	double error = 0, accurancy = eps + 1; //Погрешность решения и точость метода соответственно
	double sup; //Переменная помощник
	int IterationsCount = 0; //Количество итераций
	SetOptW();

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		accurancy = UpRelaxMethod();
		IterationsCount++;
	}

	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - u[j][i]);
			if (sup > error)
				error = sup;
		}

	sup = NevyazkaEvkl();
	TVector<double> result(4);
	result[0] = accurancy; //Точность метода
	result[1] = IterationsCount; //Количество итераций
	result[2] = error; //Погрешность решения
	result[3] = sup; //Евклидова норма невязки
	
	return result;
}


double DirichleT::NevyazkaInf()
{
	double res = 0;
	double A, x, y = yBorder[0] + k, hE, kE;
	double temp;
	hE = 1 / pow(h, 2);
	kE = 1 / pow(k, 2);
	A = -2 * (hE + kE);
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			temp = A * V[j][i] + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]) + F_Fuction(x, y);
			if (fabs(temp) > res)
				res = fabs(temp);
			x += h;
		}
		y += k;
	}
	return res;
}


double DirichleT::NevyazkaEvkl()
{
	double res = 0;
	double A, x, y = yBorder[0] + k, hE, kE;
	double temp;
	hE = 1 / pow(h, 2);
	kE = 1 / pow(k, 2);
	A = -2 * (hE + kE);
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			temp = A * V[j][i] + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]) + F_Fuction(x, y);
			res += pow(temp, 2);
			x += h;
		}
		y += k;
	}
	res = sqrt(res);
	return res;
}


void DirichleT::ChangeGrid(int N, int M)
{
	*this = DirichleT(N, M, xBorder, yBorder, TaskNumber);
}


void DirichleT::ChangeBorders(TVector<double> XBorder, TVector<double> YBorder)
{
	*this = DirichleT(n, m, XBorder, YBorder, TaskNumber);
}


void DirichleT::ChangeTask(int TASKNumber)
{
	*this = DirichleT(n, m, xBorder, yBorder, TASKNumber);
}


void DirichleT::Reset()
{
	*this = DirichleT(n, m, xBorder, yBorder, TaskNumber);
}

