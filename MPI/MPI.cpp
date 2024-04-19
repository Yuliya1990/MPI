#include <iostream>
#include <cmath>
#include <mpi.h>

double calculateTerm(int term, double x) {
    return pow(x, term) / tgamma(term + 1);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const double x = 2.0; // Вхідне значення для обчислення експоненти
    const int n = 1000; // Кількість членів ряду Тейлора

    double local_sum = 0.0;

    int next_rank = (rank + 1) % size; // Визначення наступного процесу в кільці
    int prev_rank = (rank - 1 + size) % size; // Визначення попереднього процесу в кільці
    std::cout << "next_rank " << next_rank << std::endl;
    std::cout << "prev_rank " << prev_rank << std::endl;

    // Кожен процес обчислює свій шматок ряду Тейлора

    for (int i = rank; i < n; i += size) {

        // Відправлення індексу наступному процесу
        //змінна, кількість, тип даних, номер отримувача, мітка повідомлення, комунікатор
        MPI_Send(&i, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
        // Приймання індексу від попереднього процесу
        int received_i;
        MPI_Recv(&received_i, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Продовження обчислень з отриманим індексом
        local_sum += calculateTerm(i, x);
    }
    std::cout << "Process " << rank << ": local_sum = " << local_sum << std::endl;

    double global_sum;
    //змінні що збираєсо
    //в яку змінну збираємо
    //скільки елементів буде зібрано
    // яка функція буде застосована до зібраних ел.
    // процес, що отримує результат
    // всі процеси
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double end_time = MPI_Wtime();
        std::cout << "The value of the exponent for x = " << x << " using the Taylor : " << global_sum << std::endl;
        std::cout << "Time: " << end_time - start_time << " sec" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
