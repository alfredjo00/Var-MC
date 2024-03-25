#include <stdlib.h>
#include <string.h>
#include "tools.h"
#include "metropolis.h"
/*****************************************************************************
 * Add your functions that you wanna test here (from, e.g., src/run.c.)
 *
 * Example:
 *
 * extern void get_x(...);
 *
 * ***************************************************************************/
#define SIZE 3
#include "test_main.h"

/* ************************************
 * Here follows the test we want
 * to run
 * ***********************************/
START_TEST(test_distance_between_vectors)
{
    // assert re and correct answer is within 0.001f //
    double v1[SIZE] = {1., 1., 1.};
    double v2[SIZE] = {0, -1., -1.};

    double dist_vectors12 = distance_between_vectors(v1, v2, SIZE);

    ck_assert_double_eq_tol(dist_vectors12, 3, 1e-6);


    double v3[SIZE] = {0, 5, 0.};
    double v4[SIZE] = {0, 0, 0};

    double dist_vectors34 = distance_between_vectors(v3, v4, SIZE);


    ck_assert_double_eq_tol(dist_vectors34, 5.0, 1e-6);

    // empty
}


START_TEST(test_norm)
{
    // assert re and correct answer is within 0.001f //
    double v1[SIZE] = {1., 1., 1.};
    double v2[SIZE] = {2., -1., -1.};
    double v3[SIZE] = {10, 0, 0};
    double v4[SIZE] = {0, 100, 0};

    double norm1 = vector_norm(v1, SIZE);
    double norm2 = vector_norm(v2, SIZE);
    double norm3 = vector_norm(v3, SIZE);
    double norm4 = vector_norm(v4, SIZE);

    ck_assert_double_eq_tol(norm1, sqrt(3.0), 1e-6);

    ck_assert_double_eq_tol(norm2, sqrt(6.0), 1e-6);

    ck_assert_double_eq_tol(norm3, 10.0, 1e-6);

    ck_assert_double_eq_tol(norm4, 100.0, 1e-6);

    // empty
}




int
main()
{
    test_setup("testing", "core");
    

    // Tests
    add_test(test_distance_between_vectors);
    add_test(test_norm);


    test_teardown();
    return 0;
}
