#include <mpi.h> // «адача 2 є16. Ћенточна€ вертикальна€ схема
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;
#define Main_Process 0

double** Create_matr(int _row, int _column) {// генераци€ матрицы
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
double* Create_vec(int _size) {// генераци€ вектора
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
double* Matr_to_vec(double** _matr, int _row, int _column) {// конвертаци€ матрицы в вектор
	if (_matr == NULL || _row < 1 || _column < 1)
		return NULL;

	double* _vec = new double[_row * _column];
	int id_vec = 0;
	for (int i = 0; i < _column; i++)
		for (int j = 0; j < _row; j++)
			_vec[id_vec++] = _matr[j][i];
	return _vec;
}
void Show_matr(double** _matr, double* _vec, double* _res, int _row, int _column) {// вывод матрицы в поток и в файл
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
	cout << endl << "Output matrix to file? (y/n) ";
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
int Equality_test(double* _seque, double* _paral,int _row) {// проверка идентичности результатов
	for (int i = 0; i < _row; i++)
		if (_seque[i] != _paral[i])
			return 0;
	return 1;
}

int main(int argc, char * argv[]) {
	int proc_num;			// число процессов
	int proc_rank;			// номер текущего процесса
	int numb_columns;		// число столбцов выделенных одному процессу
	int size_part_vec;		// объем данных мартицы выделенных одному процессу

	double** matr = NULL;								// матрица
	int row = 0, column = 0;							// размеры матрицы
	double* vec = NULL;									// вектор
	double* matr_as_vec = NULL;							// матрица преобразованна€ в вектор

	double* res_seque = NULL;							// итоговый вектор последовательной версии алгоритма
	double* res_paral = NULL;							// итоговый вектор параллельной версии алгоритма
	double time_seque = 0, time_paral = 0;				// врем€ работы последовательной и параллельной версии алгоритма
	double time_start_seque = 0, time_start_paral = 0;	// врем€ начала работы -//-
	double time_end_seque = 0, time_end_paral = 0;		// врем€ конца работы -//-

	double* part_matr_as_vec = NULL;					// часть матрицы принимаема€ каждым процессом
	double* part_vec = NULL;							// часть вектора -//-
	double* part_res = NULL;							// часть результата

	int* displs_m = NULL;								// массив смещений в матрице
	int* sendcounts_m = NULL;							// массив длин подматриц
	int* displs_v = NULL;								// массив смещений в векторе
	int* sendcounts_v = NULL;							// массив длин подвекторов

	int err_code = 0;
	MPI_Status stat;							// структура атрибутов сообщений дл€ "общени€" процессов
	err_code = MPI_Init(&argc, &argv);			// передача всем процессам аргументов командной строки
	if (err_code != MPI_SUCCESS) {				// проверка на успешность инициализации процессов
		cout << "Error code MPI_Init: " << err_code;
		return -1;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);	// функци€ определени€ числа процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);	// функци€ определени€ номера текущего процесса

	if (proc_rank == Main_Process) {			// код главного процесса

		cout << endl << "Number process: " << proc_num << endl;
		cout << "Enter the number of row: ";
		cin >> row;
		cout << "Enter the number of columns: ";
		cin >> column;

		res_seque = new double[row];
		res_paral = new double[row];
		for (int i = 0; i < row; i++) {
			res_seque[i] = 0;
			res_paral[i] = 0;
		}
		// генераци€ матрицы, вектора, а также представление матрицы в виде вектора
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
		// подсчет объема работы каждого процесса 
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
	}
	// Ќј„јЋќ ѕј–јЋЋ≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј

	// отправка каждому процессу массива данных и их кол-во
	MPI_Bcast(&row, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);
	MPI_Bcast(&column, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);

	numb_columns = column / proc_num;
	if (proc_rank < column % proc_num)
		numb_columns++;
	size_part_vec = numb_columns * row;

	part_matr_as_vec = new double[size_part_vec];
	part_vec = new double[numb_columns];
	part_res = new double[row];
	for (int i = 0; i < row; i++)
		part_res[i] = 0;
	
	MPI_Scatterv(matr_as_vec, sendcounts_m, displs_m, MPI_DOUBLE, part_matr_as_vec, size_part_vec, MPI_DOUBLE, Main_Process,MPI_COMM_WORLD);
	MPI_Scatterv(vec, sendcounts_v, displs_v, MPI_DOUBLE, part_vec, numb_columns, MPI_DOUBLE, Main_Process, MPI_COMM_WORLD);

	if (proc_rank == Main_Process)
		time_start_paral = MPI_Wtime();
	
	for (int i = 0; i < row; i++)
		for (int j = 0; j < numb_columns; j++)
			part_res[i] += (part_matr_as_vec[row * j + i] * part_vec[j]);

	// сложение частичных результатов в главный вектор результата
	MPI_Reduce(part_res, res_paral, row, MPI_DOUBLE, MPI_SUM, Main_Process, MPI_COMM_WORLD);
	//  ќЌ≈÷ ѕј–јЋЋ≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј

	if (proc_rank == Main_Process) {
		time_end_paral = MPI_Wtime();
		time_paral = 1000 * (time_end_paral - time_start_paral); // посчет времени работы алгоритма в миллисекудах
		cout << endl << "Spend time algorithm (Parallel version): " << time_paral << "ms" << endl << endl;

		// Ќј„јЋќ ѕќ—Ћ≈ƒќ¬ј“≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј
		time_start_seque = MPI_Wtime();
		for (int i = 0; i < row; i++)
			for (int j = 0; j < column; j++)
				res_seque[i] += (matr[i][j] * vec[j]);
		time_end_seque = MPI_Wtime();
		
		time_seque = 1000 * (time_end_seque - time_start_seque); // посчет времени работы алгоритма в миллисекудах
		cout << endl << "Spend time algorithm (Sequence version): " << time_seque << "ms" << endl;
		//  ќЌ≈÷ ѕќ—Ћ≈ƒќ¬ј“≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј

		
		// вывод результатов работы алгоритмов
		if (time_seque < time_paral)
			cout << "Sequence version faster parallel" << endl;
		else
			cout << "Parallel version faster sequence" << endl;

		if (Equality_test(res_seque, res_paral, row))
			cout << "Result sequence and parallel version identical" << endl;
		else
			cout << "Result sequence and parallel version not identical" << endl;

		// вывод в поток матрицы небольших размеров + вывод мартицы в файл
		Show_matr(matr, vec, res_paral, row, column);
		// очистка пам€ти 
		for (int i = 0; i < row; i++)
			delete[] matr[i];
		delete[] matr;
		delete[] vec;
		delete[] matr_as_vec;
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

	MPI_Finalize();// уничтощение всех MPI процессов и их св€зей
	return 0;
}