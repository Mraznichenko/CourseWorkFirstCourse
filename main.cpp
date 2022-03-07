#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void restart();
void Gauss()
{int n, i, j, k;
float d, s;
setlocale(LC_CTYPE, "ukr");
cout << "Порядок:" << endl;
cin >> n;
float **a = new float *[n];
for (i = 0; i <= n; i++)
a   [i] = new float [n];
    float **a1 = new float *[n];
for (i = 0; i <= n; i++)
    a1[i] = new float [n];
    float *b = new float [n];
    float *x = new float [n];
cout << "Введіть коефіціенти і вільні члени " << endl;
for (i = 1; i <= n; i++)
{
for (j = 1; j <= n; j++)
{
cout << "a[" << i << "," << j << "]= ";
cin >> a[i][j];
a1[i][j] = a[i][j];
}
cout << "b[" << i <<"]= ";
cin >> b[i];
}
cout << "Ваша матриця А: " << endl << endl;
for (i = 1; i <= n; i++)
{
    for (j = 1; j <= n; j++)
        cout << a[i][j] << " ";
    cout << endl;
}
for (k = 1; k <= n; k++)
{
    for (j = k + 1; j <= n; j++)
    {
        d = a[j][k] / a[k][k];
        for (i = k; i <= n; i++)
        {
        a[j][i] = a[j][i] - d * a[k][i];
        }
        b[j] = b[j] - d * b[k];
    }
}
for (k = n; k >= 1; k--)
{
    d = 0;
    for (j = k + 1; j <= n; j++)
    {
        s = a[k][j] * x[j];
        d = d + s;
    }
    x[k] = (b[k] - d) / a[k][k];
}
cout << "Корені системи: " << endl;
for( i = 1; i <= n; i++)
cout << "x[" << i << "]=" << x[i] << " " << endl;
restart();
}
void Cramer()
{
    float determinant(float matrix[3][3]);
    float determinantX1(float coefMatrix[3][3], float constTermsMatrix[3]);
    float determinantX2(float coefMatrix[3][3], float constTermsMatrix[3]);
    float determinantX3(float coefMatrix[3][3], float constTermsMatrix[3]);
    int i, j;
    float coefficientsMatrix3x3[3][3];
    float constantTermsMatrix3x1[3];
    setlocale(LC_CTYPE, "ukr");
    cout << "Введіть коефіціенти і вільні члени " << endl;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            cout << "a[" << i << "," << j << "]= ";
            cin >> coefficientsMatrix3x3[i][j];
        }
        cout << "b[" << i << "]= ";
        cin >> constantTermsMatrix3x1[i];
    }
    cout << "Ваша матриця А: " << endl << endl;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
            cout << coefficientsMatrix3x3[i][j] << " ";
        cout << endl;
    }
    float det = determinant(coefficientsMatrix3x3);
    float detX1 = determinantX1(coefficientsMatrix3x3, constantTermsMatrix3x1);
    float detX2 = determinantX2(coefficientsMatrix3x3, constantTermsMatrix3x1);
    float detX3 = determinantX3(coefficientsMatrix3x3, constantTermsMatrix3x1);

    if (det != 0)
    {
        cout << "X1 = " << (float)detX1/(float)det << endl;
        cout << "X2 = " << (float)detX2/(float)det << endl;
        cout << "X3 = " << (float)detX3/(float)det << endl;
    }
    else
        cout << "Система не має розв’язків " << endl << endl;
    restart();
}
float determinant(float matrix[3][3])
{
    float a11 = matrix[0][0];
    float a12 = matrix[0][1];
    float a13 = matrix[0][2];
    float a21 = matrix[1][0];
    float a22 = matrix[1][1];
    float a23 = matrix[1][2];
    float a31 = matrix[2][0];
    float a32 = matrix[2][1];
    float a33 = matrix[2][2];

    return (a11 * a22 * a33) + (a12 * a23 * a31) + (a13 * a21 * a32) -
           (a13 * a22 * a31) - (a11 * a23 * a32) - (a12 * a21 * a33);
}
float determinantX1(float coefMatrix[3][3], float constTermsMatrix[3])
{
    float a12 = coefMatrix[0][1];
    float a13 = coefMatrix[0][2];
    float a22 = coefMatrix[1][1];
    float a23 = coefMatrix[1][2];
    float a32 = coefMatrix[2][1];
    float a33 = coefMatrix[2][2];
    float c1 = constTermsMatrix[0];
    float c2 = constTermsMatrix[1];
    float c3 = constTermsMatrix[2];

    return (c1 * a22 * a33) + (a12 * a23 * c3) + (a13 * c2 * a32) -
           (a13 * a22 * c3) - (c1 * a23 * a32) - (a12 * c2 * a33);
}
float determinantX2(float coefMatrix[3][3], float constTermsMatrix[3])
{
    float a11 = coefMatrix[0][0];
    float a13 = coefMatrix[0][2];
    float a21 = coefMatrix[1][0];
    float a23 = coefMatrix[1][2];
    float a31 = coefMatrix[2][0];
    float a33 = coefMatrix[2][2];
    float c1 = constTermsMatrix[0];
    float c2 = constTermsMatrix[1];
    float c3 = constTermsMatrix[2];

    return (a11 * c2 * a33) + (c1 * a23 * a31) + (a13 * a21 * c3) -
           (a13 * c2 * a31) - (a11 * a23 * c3) - (c1 * a21 * a33);
}
float determinantX3(float coefMatrix[3][3], float constTermsMatrix[3])
{
    float a11 = coefMatrix[0][0];
    float a12 = coefMatrix[0][1];
    float a21 = coefMatrix[1][0];
    float a22 = coefMatrix[1][1];
    float a31 = coefMatrix[2][0];
    float a32 = coefMatrix[2][1];
    float c1 = constTermsMatrix[0];
    float c2 = constTermsMatrix[1];
    float c3 = constTermsMatrix[2];

    return (a11 * a22 * c3) + (a12 * c2 * a31) + (c1 * a21 * a32) -
           (c1 * a22 * a31) - (a11 * c2 * a32) - (a12 * a21 * c3);
}
bool converge(double xk[10], double xkp[10], int n, double eps)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
        norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
    return (sqrt(norm) < eps);
}
double okr(double x, double eps)
{
    int i = 0;
    double neweps = eps;
    while (neweps < 1)
    {
        i++;
        neweps *= 10;
    }
    int okr = pow(double(10), i);
    x = int(x * okr + 0.5) / double(okr);

    return x;
}
bool diagonal(double a[10][10], int n)
{
    int i, j, k = 1;
    double sum;
    for (i = 0; i < n; i++) {
        sum = 0;
        for (j = 0; j < n; j++)
            sum += abs(a[i][j]);
        sum -= abs(a[i][i]);
        if (sum > a[i][i])
        {
            k = 0;
            cout << a[i][i] << " < " << sum << endl;
        }
        else
        {
            cout << a[i][i] << " > " << sum << endl;
        }


    }
    return (k == 1);
}
void zeydel()
{
    setlocale(LC_ALL, "");
    double eps, a[10][10], b[10], x[10], p[10];
    int n, i, j, m = 0;
    int method;
    cout << "Введіть порядок: ";
    cin >> n;
    cout << "Оберіть точність: ";
    cin >> eps;
    cout << "Заповніть матрицю А: " << endl << endl;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            cout << "A[" << i << "][" << j << "] = ";
            cin >> a[i][j];
        }
    cout << endl << endl;
    cout << "Ваша матриця А: " << endl << endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            cout << a[i][j] << " ";
        cout << endl;
    }

    cout << endl;

    cout << "Введіть вільні члени: " << endl << endl;
    for (i = 0; i < n; i++)
    {
        cout << "В[" << i + 1 << "] = ";
        cin >> b[i];
    }
    cout << endl << endl;
    for (int i = 0; i < n; i++)
        x[i] = 1;
    if (diagonal(a, n)) {
        do
        {
            for (int i = 0; i < n; i++)
                p[i] = x[i];
            for (int i = 0; i < n; i++)
            {
                double var = 0;
                for (int j = 0; j < n; j++)
                    if(j!=i) var += (a[i][j] * x[j]);

                x[i] = (b[i] - var) / a[i][i];
            }
            m++;
        }
        while (!converge(x, p, n, eps));
        cout << "Корені системи:" << endl << endl;
        for (i = 0; i < n; i++) cout << "x" << i << " = " << okr(x[i], eps) << "" << endl;
        cout << "Ітерацій: " << m << endl;
    }
    else {
        cout << "Не выполняется умова діагоналей" << endl;
    }
    restart();
}
void del(double **A, double *X, double *Xk,   double *b, int n);
void jacobi()
{
    int n;
    int i,j;
    double **A;
    double *X;
    double *Xk;
    double *b;
    double sumd = 0;
    double eps;
    printf("Введіть порядок: \n");
    cin>>n;
    A  = new double* [n];
    X  =  new double [n];
    Xk =  new double [n];
    b  =  new double [n];
    for(i = 0; i < n; i++)
        A[i] = new double [n];
    printf("Введіть елементи системи для системи: \n");
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            scanf("%lf",&A[i][j]);
    printf("Введіть вільні члени: \n");
    for(i = 0; i < n; i++)
        scanf("%lf",&b[i]);
    printf("Задайте точність: \n");
    scanf(" %lf",&eps);
    for(i = 0; i < n; i++)
        if (A[i][i] == 0) {del( A, X, Xk, b, n); printf("Головна діагональ обнулена \n"); restart();}
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            if (i != j)
                sumd += fabs(A[i][j]);
        if (fabs(A[i][i]) < sumd)
        {
            del( A, X, Xk, b, n);
        printf("Даний метод не підходить для цієї системи \n");
        restart();
        }
        sumd = 0;
    }
    for(i = 0; i < n; i++)
    X[i]=0.0;
    int  count = 0, flag = 1;
    double x, zh = 0;
    for(i = 0; i < n; i++)
        Xk[i] = X[i];
    do
    {
        count++;
        for(i = 0; i < n; i++)
        {
            x=0;
            for(j = 0; j < n; j++)
            {
                if (i != j)
                    x += Xk[j] * A[i][j];
                if (i == j)
                    zh = A[i][j];
            }
            x = (b[i] - x) / zh;
            X[i] = x;
            if ((fabs(X[i]-Xk[i]))<=eps) flag=0;
        }
        for(i = 0; i < n; i++)
            Xk[i] = X[i];
    }while(flag);
    printf("Корені системи\n");
    for(i = 0; i < n; i++)
        printf("x[%d]=%0.5lf\n",i+1,X[i]);
    printf("iteration=%d\n",count);
    del( A, X, Xk, b, n);
    restart();
}
void del(double **A, double *X, double *Xk,   double *b, int n)
{
    for(int i = 0; i < n; i++)
        delete[]A[i];
    delete[]A;
    delete[]X;
    delete[]Xk;
    delete[]b;
}
void menu()
{
    int r;
    cout<<"Виберіть потрібний вам метод" <<endl<<"1 - Метод Гаусса"<<endl<<"2 - Метод Крамера"<<endl<<"3 - Метод Зейделя"<<endl<<"4 - Метод Якобі"<<endl;
    cin>>r;
    if (r==1)
        Gauss();
    else if(r==2)
        Cramer();
    else if(r==3)
        zeydel();
    else if(r==4)
        jacobi();
    else
        cout<<"Неприпустимий варіант"<<endl;
}
void restart()
{
    int restart;
    cout<<"Вам потрібно розв’язати ще одну СЛАР?"<<endl<<"1 - так"<<endl<<"0 - ні"<<endl;
    cin>>restart;
    if(restart)
        menu();
    else
        exit(0);
}
int main()
{
    menu();
    return 0;
}
