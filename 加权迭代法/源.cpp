#include<iostream>

using namespace std;
#define N 5  //最大迭代次数
#define ar 0.00016	//权因子

/*输入函数，输入系数矩阵A，右端B值，及初值x0*/
void input(double ** &a, double * &b, int n)
{
	int i;

	cout << "请输入系数矩阵A，以空格分隔:" << endl;
	/*输入矩阵A[][]*/
	for (i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cin >> a[i][j];
	}
	/*输入矩阵B[]*/
	cout << "请输入矩阵B，以空格分隔:" << endl;
	for (i = 0; i < n; i++)
		cin >> b[i];
}

/*利用无穷范数判断误差限*/
double max(double *x, int n)
{
	double t = fabs(x[0]);

	for (int i = 1; i < n; i++)
	{
		if (fabs(x[i] > t))
			t = fabs(x[i]);
	}
	return t;
}

/*系数矩阵A的LU分解*/
void LU_decom(double **a, double **L, double **u, int n)
{

	//max：消元时，判断矩阵是否奇异
	//sum:求解LU的列、行时使用

	int i, j, k;
	double max, sum;

	/*权因子矩阵*/
	double **I = NULL;
	I = new double *[n];

	for (i = 0; i < n; i++)
		I[i] = new double[n];

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j)I[i][j] = 1;
			else  I[i][j] = 0;
		}
	}

	/*求解LU矩阵*/
	for (k = 0; k <n; k++)
	{
		max = a[k][k];

		//LU分解
		if (fabs(max) < 1e-8) cout << "这是奇异矩阵" << endl;//此处应该跳出LU 分解
		else
		{
			//求解U的k行，L的k列
			for (i = 0; i < n; i++)
			{
				sum = 0;
				for (j = 0; j < k; j++)
					sum += L[k][j] * u[j][i];
				u[k][i] = a[k][i] + ar*I[k][j] - sum;
			}
			for (i = 0; i < n; i++)
			{
				sum = 0;
				for (j = 0; j < k; j++)
					sum += L[i][j] * u[j][k];
				L[i][k] = (a[i][k] + ar*I[i][k] - sum) / u[k][k];
			}
		}
	}

}

/*利用LU分解法求解矩阵*/
void jacobi(double **L, double **u, double *b, double *d, int n)
{
	int i, j, k;
	double sum;

	//回代，求解Y
	for (k = 0; k < n; k++)
	{

		sum = 0;
		d[k] = b[k];
		for (i = 0; i < k; i++)
			sum += L[k][i] * d[i];
		d[k] = b[k] - sum;
	}

	//回代，求解d
	for (k = n - 1; k >= 0; k--) //k>0
	{
		if (fabs(u[k][k]) < 1e-6)cout << "这是奇异矩阵" << endl;
		else
		{
			sum = 0;
			for (i = k + 1; i < n; i++)
				sum += u[k][i] * d[i];
			d[k] = (d[k] - sum) / u[k][k];
		}
	}

}

/*函数整合，总体迭代过程*/
void qiujie(int n)
{
	/*A:系数矩阵，b:右边系数，L;单位下三角矩阵，u:上三角矩阵，
	*d:jacobi法求解后的值矩阵，r:差值矩阵，x:现值
	*/
	int i, j, k;
	double sum;
	double ** a = NULL; double **L = NULL; double **u = NULL;

	double * b = NULL;
	double *x = NULL;
	double *r = NULL;
	double *d = NULL;

	a = new double *[n];            //为a开辟一个* n行空间
	L = new double *[n];			//为L开辟一个* n行空间
	u = new double *[n];		 //为u开辟一个* n行空间
	b = new double[n];				//为b开辟一个n行空间

	x = new double[n];				//为x开辟一个n行空间
	r = new double[n];				//为r开辟一个n行空间
	d = new double[n];				//为d开辟一个n行空间

	for (i = 0; i < n; i++)
	{
		a[i] = new double[n];		//为每一个*a开辟一个n行空间
		L[i] = new double[n];		//为每一个*L开辟一个n行空间
		u[i] = new double[n];		//为每一个*L开辟一个n行空间
	}

	/*调用输入函数*/
	input(a, b, n);

	/*系数矩阵A的LU分解*/
	LU_decom(a, L, u, n);
	cout << "L为：" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << L[i][j] << " ";
		cout << endl;
	}
	cout << "u为：";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << u[i][j] << " ";
		cout << endl;
	}

	/*jacobi求解值*/
	jacobi(L, u, b, d, n);

	/*记录第一次的值*/
	for (i = 0; i < n; i++)
	{
		x[i] = d[i];
	}

	/*求差值r[],*/
	for (k = 0; k < N; k++)
	{
		for (i = 0; i < n; i++)
		{
			sum = 0.0;
			for (j = 0; j < n; j++)
			{
				sum += a[i][j] * x[j];
			}
			r[i] = b[i] - sum;
		}

		/*用差值r[]的迭代，求出最后解*/
		jacobi(L, u, r, d, n);

		for (i = 0; i < n; i++)
		{
			x[i] = x[i] + d[i];
		}

		/*判断误差，并输出最后的解*/
		if (fabs(max(d, n) / max(x, n)) < 1e-5)
		{
			cout << "迭代次数为：" << k << endl;
			cout << "达到精确解,且解为：" << endl;
			for (i = 0; i < n; i++)
			{
				cout << "x[" << i + 1 << "]=" << x[i] << " , ";
			}
			cout << endl;
			break;
		}

	}
	if (k == N)
		cout << "达到最大迭代次数" << endl;
}

int main()
{
	cout << "加权迭代法求解线性方程组" << endl;
	int n;
	cout << "请输入数组的维数：" << endl;
	cin >> n;
	qiujie(n);		//调用求解函数

	system("pause");
	return 0;
}
