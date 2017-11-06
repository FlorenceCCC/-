#include<iostream>

using namespace std;
#define N 5  //����������
#define ar 0.00016	//Ȩ����

/*���뺯��������ϵ������A���Ҷ�Bֵ������ֵx0*/
void input(double ** &a, double * &b, int n)
{
	int i;

	cout << "������ϵ������A���Կո�ָ�:" << endl;
	/*�������A[][]*/
	for (i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cin >> a[i][j];
	}
	/*�������B[]*/
	cout << "���������B���Կո�ָ�:" << endl;
	for (i = 0; i < n; i++)
		cin >> b[i];
}

/*����������ж������*/
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

/*ϵ������A��LU�ֽ�*/
void LU_decom(double **a, double **L, double **u, int n)
{

	//max����Ԫʱ���жϾ����Ƿ�����
	//sum:���LU���С���ʱʹ��

	int i, j, k;
	double max, sum;

	/*Ȩ���Ӿ���*/
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

	/*���LU����*/
	for (k = 0; k <n; k++)
	{
		max = a[k][k];

		//LU�ֽ�
		if (fabs(max) < 1e-8) cout << "�����������" << endl;//�˴�Ӧ������LU �ֽ�
		else
		{
			//���U��k�У�L��k��
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

/*����LU�ֽⷨ������*/
void jacobi(double **L, double **u, double *b, double *d, int n)
{
	int i, j, k;
	double sum;

	//�ش������Y
	for (k = 0; k < n; k++)
	{

		sum = 0;
		d[k] = b[k];
		for (i = 0; i < k; i++)
			sum += L[k][i] * d[i];
		d[k] = b[k] - sum;
	}

	//�ش������d
	for (k = n - 1; k >= 0; k--) //k>0
	{
		if (fabs(u[k][k]) < 1e-6)cout << "�����������" << endl;
		else
		{
			sum = 0;
			for (i = k + 1; i < n; i++)
				sum += u[k][i] * d[i];
			d[k] = (d[k] - sum) / u[k][k];
		}
	}

}

/*�������ϣ������������*/
void qiujie(int n)
{
	/*A:ϵ������b:�ұ�ϵ����L;��λ�����Ǿ���u:�����Ǿ���
	*d:jacobi�������ֵ����r:��ֵ����x:��ֵ
	*/
	int i, j, k;
	double sum;
	double ** a = NULL; double **L = NULL; double **u = NULL;

	double * b = NULL;
	double *x = NULL;
	double *r = NULL;
	double *d = NULL;

	a = new double *[n];            //Ϊa����һ��* n�пռ�
	L = new double *[n];			//ΪL����һ��* n�пռ�
	u = new double *[n];		 //Ϊu����һ��* n�пռ�
	b = new double[n];				//Ϊb����һ��n�пռ�

	x = new double[n];				//Ϊx����һ��n�пռ�
	r = new double[n];				//Ϊr����һ��n�пռ�
	d = new double[n];				//Ϊd����һ��n�пռ�

	for (i = 0; i < n; i++)
	{
		a[i] = new double[n];		//Ϊÿһ��*a����һ��n�пռ�
		L[i] = new double[n];		//Ϊÿһ��*L����һ��n�пռ�
		u[i] = new double[n];		//Ϊÿһ��*L����һ��n�пռ�
	}

	/*�������뺯��*/
	input(a, b, n);

	/*ϵ������A��LU�ֽ�*/
	LU_decom(a, L, u, n);
	cout << "LΪ��" << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << L[i][j] << " ";
		cout << endl;
	}
	cout << "uΪ��";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << u[i][j] << " ";
		cout << endl;
	}

	/*jacobi���ֵ*/
	jacobi(L, u, b, d, n);

	/*��¼��һ�ε�ֵ*/
	for (i = 0; i < n; i++)
	{
		x[i] = d[i];
	}

	/*���ֵr[],*/
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

		/*�ò�ֵr[]�ĵ������������*/
		jacobi(L, u, r, d, n);

		for (i = 0; i < n; i++)
		{
			x[i] = x[i] + d[i];
		}

		/*�ж�����������Ľ�*/
		if (fabs(max(d, n) / max(x, n)) < 1e-5)
		{
			cout << "��������Ϊ��" << k << endl;
			cout << "�ﵽ��ȷ��,�ҽ�Ϊ��" << endl;
			for (i = 0; i < n; i++)
			{
				cout << "x[" << i + 1 << "]=" << x[i] << " , ";
			}
			cout << endl;
			break;
		}

	}
	if (k == N)
		cout << "�ﵽ����������" << endl;
}

int main()
{
	cout << "��Ȩ������������Է�����" << endl;
	int n;
	cout << "�����������ά����" << endl;
	cin >> n;
	qiujie(n);		//������⺯��

	system("pause");
	return 0;
}
