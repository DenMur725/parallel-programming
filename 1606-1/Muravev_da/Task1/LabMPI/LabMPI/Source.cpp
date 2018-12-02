#include <mpi.h> // «адача 1 є10 —умма элементов матрицы
#include <ctime>
#include <iostream>
#include <fstream>
using namespace std;
#define Main_Process 0

double* Create_matr_as_vec(int _size) {// генераци€ матрицы в векторе
	if ( _size < 1)
		return NULL;
	int a = 0, b = 0;
	cout << endl << "Enter range [a,b] for matrix" << endl;
	cin >> a >> b;
	double*  _matr_as_vec = new double[_size];
	for (int i = 0; i < _size; i++)
			_matr_as_vec[i] = rand() % (b - a + 1) + a;
	return  _matr_as_vec;
}
void Show_matr(double* _matr_as_vec, int _row, int _column) {// вывод матрицы в поток и в файл
	if (_matr_as_vec == NULL || _row < 1 || _column < 1)
		return;
	int delta = 7;
	int id_vec = 0;
	if (_row < 18 || _column < 8) {
		cout << "Matrix [" << _row << "x" << _column << "]" << endl;

		for (int i = 0; i < _row; i++) {
			for (int j = 0; j < _column; j++) {
				cout.width(delta);
				cout << _matr_as_vec[id_vec++];
			}
			cout << endl << endl;
		}
	}
	id_vec = 0;
	cout << endl << "Output matrix to file? (y/n) ";
	char answer = 'n';
	cin >> answer;
	if (answer == 'y') {
		ofstream strm("Matrix.txt");
		strm << "Matrix [" << _row << "x" << _column << "]" << endl;
		for (int i = 0; i < _row; i++) {
			for (int j = 0; j < _column; j++) {
				strm.width(delta);
				strm << _matr_as_vec[id_vec++];
			}
			strm << endl << endl;
		}
		strm.close();
	}
}
int Equality_test(double _seque, double _paral) {// проверка идентичности результатов
	if (_seque != _paral)
		return 0;
	return 1;
}

int main(int argc, char * argv[]) {
	int proc_num;			// число процессов
	int proc_rank;			// номер текущего процесса
	int size_part_vec;		// объем данных мартицы выделенных одному процессу
	
	double* matr_as_vec = NULL;							// матрица преобразованна€ в вектор
	int row = 0, column = 0;							// размеры матрицы
	int size = 0;										// длина вектора

	double sum_seque = 0;								// итогова€ сумма последовательной версии алгоритма
	double sum_paral = 0;								// итогова€ сумма параллельной версии алгоритма

	double time_seque = 0, time_paral = 0;				// врем€ работы последовательной и параллельной версии алгоритма
	double time_start_seque = 0, time_start_paral = 0;	// врем€ начала работы -//-
	double time_end_seque = 0, time_end_paral = 0;		// врем€ конца работы -//-

	double* part_matr_as_vec = NULL;					// часть матрицы принимаема€ каждым процессом
	double part_sum = 0;								// частична€ сумма

	int* displs = NULL;									// массив смещений в векторе
	int* sendcounts = NULL;								// массив длин подвекторов

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
		size = row * column;

		// генераци€ матрицы, вектора, а также представление матрицы в виде вектора
		srand((unsigned)time(NULL));
		matr_as_vec = Create_matr_as_vec(size);
		if (matr_as_vec == NULL) {
			cout << "Error! Incorrect input data for matrix";
			MPI_Finalize();
			return -1;
		}
		// подсчет объема работы каждого процесса 
		size_part_vec = size / proc_num;
		int ost = size % proc_num;

		displs = new int[proc_num];
		sendcounts = new int[proc_num];
		for (int i = 0; i < ost; i++)
			sendcounts[i] = size_part_vec + 1;
		for (int i = ost; i < proc_num; i++)
			sendcounts[i] = size_part_vec;
		displs[0] = 0;
		for (int i = 1; i < proc_num; i++) {
			displs[i] = displs[i - 1] + sendcounts[i - 1];
		}
	}
	// Ќј„јЋќ ѕј–јЋЋ≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј

	// отправка каждому процессу массива данных и их кол-во
	MPI_Bcast(&size, 1, MPI_INT, Main_Process, MPI_COMM_WORLD);

	size_part_vec = size / proc_num;
	if (proc_rank < size % proc_num)
		size_part_vec++;
	part_matr_as_vec = new double[size_part_vec];

	MPI_Scatterv(matr_as_vec, sendcounts, displs, MPI_DOUBLE, part_matr_as_vec, size_part_vec, MPI_DOUBLE, Main_Process, MPI_COMM_WORLD);

	if (proc_rank == Main_Process)
		time_start_paral = MPI_Wtime();

	for (int i = 0; i < size_part_vec; i++)
		part_sum += part_matr_as_vec[i];

	// сложение частичных результатов в итоговую сумму
	MPI_Reduce(&part_sum, &sum_paral, 1, MPI_DOUBLE, MPI_SUM, Main_Process, MPI_COMM_WORLD);

	//  ќЌ≈÷ ѕј–јЋЋ≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј

	if (proc_rank == Main_Process) {
		time_end_paral = MPI_Wtime();
		time_paral = 1000 * (time_end_paral - time_start_paral); // посчет времени работы алгоритма в миллисекудах
		cout << endl << "Spend time algorithm (Parallel version): " << time_paral << "ms" << endl << endl;

		// Ќј„јЋќ ѕќ—Ћ≈ƒќ¬ј“≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј
		time_start_seque = MPI_Wtime();
		for (int i = 0; i < size; i++)
			sum_seque += matr_as_vec[i];
		time_end_seque = MPI_Wtime();

		time_seque = 1000 * (time_end_seque - time_start_seque); // посчет времени работы алгоритма в миллисекудах
		cout << endl << "Spend time algorithm (Sequence version): " << time_seque << "ms" << endl;
		//  ќЌ≈÷ ѕќ—Ћ≈ƒќ¬ј“≈Ћ№Ќќ√ќ јЋ√ќ–»“ћј


		// вывод результатов работы алгоритмов
		if (time_seque < time_paral)
			cout << "Sequence version faster parallel" << endl;
		else
			cout << "Parallel version faster sequence" << endl;

		if (Equality_test(sum_seque, sum_paral))
			cout << "Result sequence and parallel version identical" << endl;
		else
			cout << "Result sequence and parallel version not identical" << endl;

		// вывод в поток матрицы небольших размеров + вывод мартицы в файл
		Show_matr(matr_as_vec, row, column);
		// очистка пам€ти 
		delete[] matr_as_vec;
		delete[] displs;
		delete[] sendcounts;
	}
	delete[] part_matr_as_vec;

	MPI_Finalize();// уничтощение всех MPI процессов и их св€зей
	return 0;
}