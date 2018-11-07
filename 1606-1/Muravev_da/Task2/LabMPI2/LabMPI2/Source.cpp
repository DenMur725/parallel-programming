#include <mpi.h> // ������ 2 �16. ��������� ������������ �����
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;
#define Main_Process		0
#define Port_Size_Vector	1
#define Port_Size_Column	2
#define Port_Matrix			4
#define Port_Vector			8
#define Port_Result			16

int proc_num;			// ����� ���������
int proc_rank;			// ����� �������� ��������
int numb_columns;		// ����� �������� ���������� ������ ��������
int size_part_vec;		// ����� ������ ������� ���������� ������ ��������


double** Create_matr(int _row, int _column) {// ��������� �������
	if (_row < 1 || _column < 1)
		return NULL;

	double** _matr = new double*[_row];
	for (int i = 0; i < _row; i++)
		_matr[i] = new double[_column];
	int a = 0, b = 0;
	cout << endl << "Enter range [a,b] for matrix" << endl;
	cin >> a >> b;
	for (int i = 0; i < _row; i++)
		for (int j = 0; j < _column; j++)
			_matr[i][j] = rand() % (b - a + 1) + a;
	return _matr;
}
double* Create_vec(int _size) {// ��������� �������
	if (_size < 1)
		return NULL;

	double* _vec = new double[_size];
	int a = 0, b = 0;
	cout << endl << "Enter range [a,b] for vector" << endl;
	cin >> a >> b;
	for (int i = 0; i < _size; i++)
		_vec[i] = rand() % (b - a + 1) + a;
	return _vec;
}
double* Matr_to_vec(double** _matr, int _row, int _column) {// ����������� ������� � ������
	if (_matr == NULL || _row < 1 || _column < 1)
		return NULL;

	double* _vec = new double[_row * _column];
	int id_vec = 0;
	for (int i = 0; i < _column; i++)
		for (int j = 0; j < _row; j++)
			_vec[id_vec++] = _matr[j][i];
	return _vec;
}
void Show_matr(double** _matr, double* _vec, double* _res, int _row, int _column) {// ����� ������� � ����� � � ����
	if (_matr == NULL || _row < 1 || _column < 1)
		return;
	int delta = 7;
	if (_row < 18 || _column < 8) {
		cout << "Matrix [" << _row << "x" << _column << "]" <<endl;

		for (int i = 0; i < _row; i++) {
			for (int j = 0; j < _column; j++) {
				cout.width(delta);
				cout << _matr[i][j];
			}
			cout << endl << endl;
		}
		cout << endl << endl << "Vector [" <<_column << "]" << endl;
		for (int i = 0; i < _column; i++) {
			cout.width(delta);
			cout << _vec[i];
		}
		cout << endl << endl << "Result [" << _row << "]" << endl;
		for (int i = 0; i < _row; i++) {
			cout.width(delta);
			cout << _res[i];
		}
	}
	cout << "Output matrix to file? (y/n) ";
	char answer = 'n';
	cin >> answer;
	if (answer == 'y') {
		ofstream strm("Matrix X Vector.txt");
		strm << "Matrix [" << _row << "x" << _column << "]" << endl;
		for (int i = 0; i < _row; i++) {
			for (int j = 0; j < _column; j++) {
				strm.width(delta);
				strm << _matr[i][j];
			}
			strm << endl << endl;
		}
		strm << endl << endl << "Vector [" << _column << "]" << endl;
		for (int i = 0; i < _column; i++) {
			strm.width(delta);
			strm << _vec[i];
		}
		strm << endl << endl << "Result [" << _row << "]" << endl;
		for (int i = 0; i < _row; i++) {
			strm.width(delta);
			strm << _res[i];
		}
		strm.close();
	}
}
int Equality_test(double* _seque, double* _paral,int _row) {// �������� ������������ �����������
	for (int i = 0; i < _row; i++)
		if (_seque[i] != _paral[i])
			return 0;
	return 1;
}

int main(int argc, char * argv[]) {
	int err_code;
	MPI_Status stat;							// ��������� ��������� ��������� ��� "�������" ���������
	err_code = MPI_Init(&argc, &argv);			// �������� ���� ��������� ���������� ��������� ������
	if (err_code != MPI_SUCCESS) {				// �������� �� ���������� ������������� ���������
		cout << "Error code MPI_Init: " << err_code;
		return -1;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);	// ������� ����������� ����� ���������
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);	// ������� ����������� ������ �������� ��������

	if (proc_num < 1) {
		cout << "Error! The number of processes must be at least 1";
		return -1;
	}
	if (proc_rank == Main_Process) {						// ��� �������� ��������

		double** matr = NULL;								// �������
		int row = 0, column = 0;							// ������� �������
		int size_matr = 0;									// ����� ��������� � �������
		double* vec = NULL;									// ������
		double* matr_as_vec = NULL;							// ������� ��������������� � ������

		double* res_seque = NULL;							// �������� ������ ���������������� ������ ���������
		double* res_paral = NULL;							// �������� ������ ������������ ������ ���������
		double time_seque = 0, time_paral = 0;				// ����� ������ ���������������� � ������������ ������ ���������
		double time_start_seque = 0, time_start_paral = 0;	// ����� ������ ������ -//-
		double time_end_seque = 0, time_end_paral = 0;		// ����� ����� ������ -//-

		cout << endl << "Number process: " << proc_num << endl;
		cout << "Enter the number of row: ";
		cin >> row;
		cout << "Enter the number of columns: ";
		cin >> column;
		size_matr = row * column;

		res_seque = new double[row];
		res_paral = new double[row];
		double* part_sum = new double[row];
		for (int i = 0; i < row; i++) {
			res_seque[i] = 0;
			res_paral[i] = 0;
		}
		// ��������� �������, �������, � ����� ������������� ������� � ���� �������
		srand((unsigned)time(NULL));
		matr = Create_matr(row, column);
		if (matr == NULL) {
			cout << "Error! Incorrect input data for matrix";
			return -1;
		}
		vec = Create_vec(column);
		if (vec == NULL) {
			cout << "Error! Incorrect input data for vector";
			return -1;
		}
		matr_as_vec = Matr_to_vec(matr, row, column);
		if (matr_as_vec == NULL) {
			cout << "Error! Incorrect input data for matrix as vector";
			return -1;
		}

		// ������ ����������������� ���������
		time_start_seque = MPI_Wtime();
		for (int i = 0; i < row; i++) 
			for (int j = 0; j < column; j++)
				res_seque[i] += (matr[i][j] * vec[j]);
		time_end_seque = MPI_Wtime();

		// ������ ������� ������ ��������� � ������������
		time_seque = 1000 * (time_end_seque - time_start_seque);
		cout << endl << "Spend time algorithm (Sequence version): " << time_seque << "ms" << endl;
		// ����� ����������������� ���������

		
		// ������� ������ ������ ������� ��������
		numb_columns = column / proc_num;
		size_part_vec = numb_columns * row;

		// ������ ������������� ���������
		// �������� ������� �������� ������� ������ � �� ���-��

		for (int i = 1; i < proc_num; i++)
			MPI_Send(&size_part_vec, 1, MPI_INT, i, Port_Size_Vector, MPI_COMM_WORLD);
		for (int i = 1; i < proc_num; i++)
			MPI_Send(&row, 1, MPI_INT, i, Port_Size_Column, MPI_COMM_WORLD);
		
		for (int i = 1; i < proc_num; i++)
			MPI_Send(matr_as_vec + size_part_vec * (i - 1), size_part_vec, MPI_DOUBLE, i, Port_Matrix, MPI_COMM_WORLD);
		for (int i = 1; i < proc_num; i++)
			MPI_Send(vec + numb_columns * (i - 1), numb_columns, MPI_DOUBLE, i, Port_Vector, MPI_COMM_WORLD);

		time_start_paral = MPI_Wtime();
		// ������� ���������� ��������  (column % proc_num)

		for (int i = 0; i < row; i++)
			for (int j = numb_columns * (proc_num - 1); j < column; j++)
				res_paral[i] += (matr_as_vec[row * j + i] * vec[j]);

		// �������� �� ���� ��������� ������������ ��������� ���� ������������
		for (int i = 1; i < proc_num; i++) {
			MPI_Recv(part_sum, row, MPI_DOUBLE, i, Port_Result, MPI_COMM_WORLD, &stat);
			for (int j = 0; j < row; j++)
				res_paral[j] += part_sum[j];
		}
		time_end_paral = MPI_Wtime();

		// ������ ������� ������ ��������� � ������������
		time_paral = 1000 * (time_end_paral - time_start_paral);
		cout << endl << "Spend time algorithm (Parallel version): " << time_paral << "ms" << endl << endl;
		// ����� ������������� ���������

		// ����� ����������� ������ ����������
		if (time_seque < time_paral)
			cout << "Sequence version faster parallel" << endl;
		else
			cout << "Parallel version faster sequence" << endl;

		if (Equality_test(res_seque, res_paral, row))
			cout << "Result sequence and parallel version identical" << endl;
		else
			cout << "Result sequence and parallel version not identical" << endl;

		// ����� � ����� ������� ��������� �������� + ����� ������� � ����
		Show_matr(matr, vec, res_paral, row, column);
		// ������� ������
		for (int i = 0; i < row; i++)
			delete[] matr[i];
		delete[] matr;
		delete[] vec;
		delete[] matr_as_vec;
		delete[] res_seque;
		delete[] res_paral;

	}
	else {
		// ����� ������ �� �������� ��������
		int _row = 0;
		double* _matr_as_vec = NULL;
		double* _vec = NULL;
		double* _res = NULL;

		MPI_Recv(&size_part_vec, 1, MPI_INT, Main_Process, Port_Size_Vector, MPI_COMM_WORLD, &stat);
		MPI_Recv(&_row, 1, MPI_INT, Main_Process, Port_Size_Column, MPI_COMM_WORLD, &stat);

		numb_columns = size_part_vec / _row;
		_matr_as_vec = new double[size_part_vec];
		_vec = new double[numb_columns];
		_res = new double[_row];

		MPI_Recv(_matr_as_vec, size_part_vec, MPI_DOUBLE, Main_Process, Port_Matrix, MPI_COMM_WORLD, &stat);
		MPI_Recv(_vec, numb_columns, MPI_DOUBLE, Main_Process, Port_Vector, MPI_COMM_WORLD, &stat);
		// ������� ��������� �����
		for (int i = 0; i < _row; i++)
			_res[i] = 0;

		for (int i = 0; i < _row; i++)
			for (int j = 0; j < numb_columns; j++)
				_res[i] += (_matr_as_vec[_row * j + i] * _vec[j]);

		// �������� � ������� ������� ��������� �����
		MPI_Send(_res, _row, MPI_DOUBLE, Main_Process, Port_Result, MPI_COMM_WORLD);
		// ������� ������
		delete[] _matr_as_vec;
		delete[] _vec;
		delete[] _res;
	}
	MPI_Finalize();// ����������� ���� MPI ��������� � �� ������
	return 0;
}