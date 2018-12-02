#include <mpi.h> // ������ 2 �16. ��������� ������������ �����
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;
#define Main_Process 0

double* Create_vec(int _size, int _a = 0, int _b = 1) {// ��������� �������
	if (_size < 1)
		return NULL;
	double* _vec = new double[_size];
	for (int i = 0; i < _size; i++)
		_vec[i] = rand() % (_b - _a + 1) + _a;
	return _vec;
}
void Show_matr(double* _matr_as_vec, double* _vec, double* _res, int _row, int _column) {// ����� ������� � �������� � ����� � � ����
	if (_matr_as_vec == NULL || _row < 1 || _column < 1)
		return;
	int id_vec = 0;
	int delta = 7;
	if (_row < 18 || _column < 8) {
		cout << "Matrix [" << _row << "x" << _column << "]" <<endl;
		
		for (int i = 0; i < _row; i++) {
			for (int j = 0; j < _column; j++) {
				cout.width(delta);
				cout << _matr_as_vec[id_vec++];
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
	id_vec = 0;
	cout << endl << "Output matrix to file? (y/n) ";
	char answer = 'n';
	cin >> answer;
	if (answer == 'y') {
		ofstream strm("Matrix X Vector.txt");
		strm << "Matrix [" << _row << "x" << _column << "]" << endl;
		for (int i = 0; i < _row; i++) {
			for (int j = 0; j < _column; j++) {
				strm.width(delta);
				strm << _matr_as_vec[id_vec++];
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
int Equality_test(double* _seque, double* _paral, int _row) {// �������� ������������ �����������
	for (int i = 0; i < _row; i++)
		if (_seque[i] != _paral[i])
			return 0;
	return 1;
}


int main(int argc, char * argv[]) {
	int proc_num;			// ����� ���������
	int proc_rank;			// ����� �������� ��������
	int numb_columns;		// ����� �������� ���������� ������ ��������
	int size_part_vec;		// ����� ������ ������� ���������� ������ ��������

	int row = 0, column = 0, size = 0;					// ������� �������
	double* vec = NULL;									// ������
	double* matr_as_vec = NULL;							// ������� ��������������� � ������
	char* buff_col = NULL;								// ����� ��������

	double* res_seque = NULL;							// �������� ������ ���������������� ������ ���������
	double* res_paral = NULL;							// �������� ������ ������������ ������ ���������
	double time_seque = 0, time_paral = 0;				// ����� ������ ���������������� � ������������ ������ ���������
	double time_start_seque = 0, time_start_paral = 0;	// ����� ������ ������ -//-
	double time_end_seque = 0, time_end_paral = 0;		// ����� ����� ������ -//-

	double* part_matr_as_vec = NULL;					// ����� ������� ����������� ������ ���������
	double* part_vec = NULL;							// ����� ������� -//-
	double* part_res = NULL;							// ����� ����������

	int* displs_m = NULL;								// ������ �������� � �������
	int* sendcounts_m = NULL;							// ������ ���� ���������
	int* displs_v = NULL;								// ������ �������� � �������
	int* sendcounts_v = NULL;							// ������ ���� �����������

	int err_code = 0;
	MPI_Status stat;							// ��������� ��������� ��������� ��� "�������" ���������
	err_code = MPI_Init(&argc, &argv);			// �������� ���� ��������� ���������� ��������� ������
	if (err_code != MPI_SUCCESS) {				// �������� �� ���������� ������������� ���������
		cout << "Error code MPI_Init: " << err_code;
		return -1;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);	// ������� ����������� ����� ���������
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);	// ������� ����������� ������ �������� ��������

	if (proc_rank == Main_Process) {			// ��� �������� ��������

		cout << endl << "Number process: " << proc_num << endl;
		cout << "Enter the number of row: ";
		cin >> row;
		cout << "Enter the number of columns: ";
		cin >> column;
		size = row * column;

		res_seque = new double[row];
		res_paral = new double[row];
		for (int i = 0; i < row; i++) {
			res_seque[i] = 0;
			res_paral[i] = 0;
		}
		// ��������� ������� � �������
		srand((unsigned)time(NULL));
		int a = 0, b = 1;
		cout << endl << "Enter range [a,b] for matrix" << endl;
		cin >> a >> b;
		matr_as_vec = Create_vec(size, a, b);
		if (matr_as_vec == NULL) {
			cout << "Error! Incorrect input data for matrix";
			MPI_Finalize();
			return -1;
		}
		cout << "Enter range [a,b] for vector" << endl;
		cin >> a >> b;
		vec = Create_vec(column, a, b);
		if (vec == NULL) {
			cout << "Error! Incorrect input data for vector";
			MPI_Finalize();
			return -1;
		}

		// ������� ������ ������ ������� �������� 
		numb_columns = column / proc_num;
		int ost = column % proc_num;

		displs_m = new int[proc_num];
		sendcounts_m = new int[proc_num];
		displs_v = new int[proc_num];
		sendcounts_v = new int[proc_num];

		for (int i = 0; i < ost; i++)
			sendcounts_v[i] = numb_columns + 1;
		for (int i = ost; i < proc_num; i++)
			sendcounts_v[i] = numb_columns;

		sendcounts_m[0] = sendcounts_v[0] * row;
		displs_v[0] = 0;
		displs_m[0] = 0;
		for (int i = 1; i < proc_num; i++) {
			sendcounts_m[i] = sendcounts_v[i] * row;
			displs_v[i] = displs_v[i - 1] + sendcounts_v[i - 1];
			displs_m[i] = displs_m[i - 1] + sendcounts_m[i - 1];
		}

		MPI_Datatype col, col_type;
		MPI_Type_vector(row, 1, column, MPI_DOUBLE, &col);
		MPI_Type_commit(&col);
		MPI_Type_create_resized(col, 0, 1 * sizeof(double), &col_type);
		MPI_Type_commit(&col_type);

		int position = 0;
		int size_buf = size * sizeof(double);
		buff_col = new char[size_buf];
		MPI_Pack(matr_as_vec, column, col_type, buff_col, size_buf, &position, MPI_COMM_WORLD);
		MPI_Type_free(&col);
		MPI_Type_free(&col_type);

	}
	// ������ ������������� ���������

	// �������� ������� �������� ������� ������ � �� ���-��
	MPI_Bcast(&row, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);
	MPI_Bcast(&column, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);

	numb_columns = column / proc_num;
	int ost = column % proc_num;
	if (proc_rank < ost)
		numb_columns++;
	size_part_vec = row * numb_columns;
	part_matr_as_vec = new double[size_part_vec];
	part_vec = new double[numb_columns];
	part_res = new double[row];
	for (int i = 0; i < row; i++)
		part_res[i] = 0;
	
	MPI_Scatterv(buff_col, sendcounts_m, displs_m, MPI_DOUBLE, part_matr_as_vec, size_part_vec, MPI_DOUBLE, Main_Process, MPI_COMM_WORLD);
	MPI_Scatterv(vec, sendcounts_v, displs_v, MPI_DOUBLE, part_vec, numb_columns, MPI_DOUBLE, Main_Process, MPI_COMM_WORLD);

	if (proc_rank == Main_Process)
		time_start_paral = MPI_Wtime();

	for (int i = 0; i < row; i++)
		for (int j = 0; j < numb_columns; j++)
			part_res[i] += (part_matr_as_vec[row * j + i] * part_vec[j]);
	
	// �������� ��������� ����������� � ������� ������ ����������
	MPI_Reduce(part_res, res_paral, row, MPI_DOUBLE, MPI_SUM, Main_Process, MPI_COMM_WORLD);

	// ����� ������������� ���������

	if (proc_rank == Main_Process) {
		time_end_paral = MPI_Wtime();
		time_paral = 1000 * (time_end_paral - time_start_paral); // ������ ������� ������ ��������� � ������������
		cout << endl << "Spend time algorithm (Parallel version): " << time_paral << "ms";

		// ������ ����������������� ���������
		time_start_seque = MPI_Wtime();
		int id_vec = 0;
		for (int i = 0; i < row; i++)
			for (int j = 0; j < column; j++)
				res_seque[i] += (matr_as_vec[id_vec++] * vec[j]);
		time_end_seque = MPI_Wtime();
		
		time_seque = 1000 * (time_end_seque - time_start_seque); // ������ ������� ������ ��������� � ������������
		cout << endl << "Spend time algorithm (Sequence version): " << time_seque << "ms" << endl;
		// ����� ����������������� ���������

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
		Show_matr(matr_as_vec, vec, res_paral, row, column);
		// ������� ������ 
		delete[] vec;
		delete[] matr_as_vec;
		delete[] buff_col;
		delete[] res_seque;
		delete[] res_paral;
		delete[] displs_m;
		delete[] sendcounts_m;
		delete[] displs_v;
		delete[] sendcounts_v;
	}
	delete[] part_matr_as_vec;
	delete[] part_vec;
	delete[] part_res;
	
	MPI_Finalize();// ����������� ���� MPI ��������� � �� ������
	return 0;
}