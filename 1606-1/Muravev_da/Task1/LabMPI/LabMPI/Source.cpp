#include <mpi.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace std;
#define Main_Process		0
#define Port_Size_Vector	2
#define Port_Vector			4
#define Port_Part_Sum		8

int proc_num;			// ����� ���������
int proc_rank;			// ����� �������� ��������
int size_vec_for_proc;	// ����� ������� ����������� ����� ��������� ���������� ��� 1 ��������
double part_sum;		// ����� ��������� ���������� ��� 1 ��������

double** Create_matr(int _size_row, int _size_column) {// ��������� �������
	if (_size_row < 1 || _size_column < 1)
		return NULL;

	srand((unsigned)time(NULL));
	double** _matr = new double*[_size_row];
	for (int i = 0; i < _size_row; i++)
		_matr[i] = new double[_size_column];
	
	for (int i = 0; i < _size_row; i++)
		for (int j = 0; j < _size_column; j++)
			_matr[i][j] = (rand() % 200)  - 99;
	return _matr;
}

double* Matr_to_vec(double** _matr, int _size_row, int _size_column) {// ����������� ������� � ������
	if (_matr == NULL || _size_row < 1 || _size_column < 1)
		return NULL;

	double* _vec = new double[_size_row * _size_column];
	int id_vec = 0;
	for (int i = 0; i < _size_row; i++)
		for (int j = 0; j < _size_column; j++)
			_vec[id_vec++] = _matr[i][j];
	return _vec;
}

void Show_matr(double** _matr, int _size_row, int _size_column) {// ����� ������� � ����� � � ����
	if (_matr == NULL || _size_row < 1 || _size_column < 1)
		return;
	if (_size_row < 21 || _size_column < 11) {
		cout << "Matrix " << _size_row << "x" << _size_column << endl;
		
		for (int i = 0; i < _size_row; i++) {
			for (int j = 0; j < _size_column; j++) {
				cout.width(7);
				cout << _matr[i][j];
			}
			cout << endl << endl;
		}
	}
	cout << "Output matrix to file? (y/n) ";
	char answer = 'n';
	cin >> answer;
	if (answer == 'y') {
		ofstream strm("Matrix.txt");
		strm << "Matrix " << _size_row << "x" << _size_column << endl;
		for (int i = 0; i < _size_row; i++) {
			for (int j = 0; j < _size_column; j++) {
				strm.width(7);
				strm << _matr[i][j];
			}
			strm << endl << endl;
		}
		strm.close();
	}
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
		printf("Error! The number of processes must be at least 1");
		return -1;
	}
	if (proc_rank == Main_Process) {						// ��� �������� ��������

		double** matr = NULL;								// �������
		int size_row = 0, size_column = 0;					// ������� �������

		double* matr_as_vec = NULL;							// ������� ��������������� � ������
		int size_vec;										// ����� ������� ����������� ��� �������� �������

		double sum_seque = 0, sum_paral = 0;				// ����� ��������� ���������������� � ������������ ������ ���������
		double time_seque = 0, time_paral = 0;				// ����� ������ -//-
		double time_start_seque = 0, time_start_paral = 0;	// ����� ������ ������ -//-
		double time_end_seque = 0, time_end_paral = 0;		// ����� ����� ������ -//-

		// ���� �������� �������
		cout << "Input size of row: ";
		cin >> size_row;
		cout << "Input size of column: ";
		cin >> size_column;

		// ��������� ������� � �������������� �� � ������
		matr = Create_matr(size_row, size_column);
		if (matr == NULL) {
			printf("Error! Incorrect input data for matrix");
			return -1;
		}
		size_vec = size_row * size_column;
		matr_as_vec = Matr_to_vec(matr, size_row, size_column);
		if (matr == NULL) {
			printf("Error! Incorrect input data for vector");
			return -1;
		}

		// ����� � ��������� ������ ������� ��������� �������� + ����� ������� � ����
		Show_matr(matr, size_row, size_column);

		// ������ ����������������� ���������
		time_start_seque = MPI_Wtime();
		for (int i = 0; i < size_row; i++)
			for (int j = 0; j < size_column; j++)
				sum_seque += matr[i][j];
		time_end_seque = MPI_Wtime();

		// ������ ������� ������ ��������� � ������������
		time_seque = 1000 * (time_end_seque - time_start_seque);
		cout << endl << "Spend time algorithm (Sequence version): " << time_seque << "ms";
		cout << endl << "Sum of all elements in matrix (Sequence version): " << sum_seque << endl;
		// ����� ����������������� ���������

		cout << endl << "Number process: " << proc_num << endl;
		// ������� ������ ������ ������� ��������
		size_vec_for_proc = size_vec / proc_num;

		// ������ ������������� ���������
		// �������� ������� �������� ������� ������ � �� ���-��
		for (int i = 1; i < proc_num; i++)
			MPI_Send(&size_vec_for_proc, 1, MPI_INT, i, Port_Size_Vector, MPI_COMM_WORLD);
		for (int i = 1; i < proc_num; i++)
			MPI_Send(matr_as_vec + size_vec_for_proc * (i - 1), size_vec_for_proc, MPI_DOUBLE, i, Port_Vector, MPI_COMM_WORLD);

		time_start_paral = MPI_Wtime();
		// ������� ���������� ��������� �������  (size_vec % proc_num)
		part_sum = 0;
		for (int i = size_vec_for_proc * (proc_num - 1); i < size_vec; i++)
			part_sum += matr_as_vec[i];
		sum_paral = part_sum;

		// �������� �� ���� ��������� ������������ ��������� ����
		for (int i = 1; i < proc_num; i++) {
			MPI_Recv(&part_sum, 1, MPI_DOUBLE, i, Port_Part_Sum, MPI_COMM_WORLD, &stat);
			sum_paral += part_sum;
		}
		time_end_paral = MPI_Wtime();

		// ������ ������� ������ ��������� � ������������
		time_paral = 1000 * (time_end_paral - time_start_paral);
		cout << endl << "Spend time algorithm (Parallel version): " << time_paral << "ms";
		cout << endl << "Sum of all elements in matrix (Parallel version): " << sum_paral << endl << endl;
		// ����� ������������� ���������

		// ����� ����������� ������ ����������
		if (time_seque < time_paral)
			cout << "Sequence version faster parallel" << endl;
		else
			cout << "Parallel version faster sequence" << endl;

		if (sum_seque == sum_paral)
			cout << "Sum sequence and parallel version identical" << endl;
		else
			cout << "Sum sequence and parallel version not identical" << endl;

		// ������� ������
		for (int i = 0; i < size_row; i++)
			delete matr[i];
		delete matr;
		delete matr_as_vec;
	}
	else { // ��� ��������������� ���������

		// ����� ������ �� �������� ��������
		MPI_Recv(&size_vec_for_proc, 1, MPI_INT, Main_Process, Port_Size_Vector, MPI_COMM_WORLD, &stat);
		double* vec = new double[size_vec_for_proc];
		MPI_Recv(vec, size_vec_for_proc, MPI_DOUBLE, Main_Process, Port_Vector, MPI_COMM_WORLD, &stat);

		// ������� ��������� �����
		part_sum = 0;
		for (int i = 0; i < size_vec_for_proc; i++)
			part_sum += vec[i];

		// �������� � ������� ������� ��������� �����
		MPI_Send(&part_sum, 1, MPI_DOUBLE, Main_Process, Port_Part_Sum, MPI_COMM_WORLD);
		delete vec;
	}
	MPI_Finalize();// ����������� ���� MPI ��������� � �� ������
	return 0;
}