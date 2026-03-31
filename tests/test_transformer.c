/*
 * Unit Tests for Transformer Models
 * Tests ideal transformer accuracy
 */

#include "../mna_solver.h"
#include "test_framework.h"
#include <stdio.h>
#include <stdlib.h>

/* Test 1: Ideal transformer turns ratio retrieval */
static int test_ideal_transformer_creation(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* Create ideal transformer with 2:1 turns ratio */
    ComponentHandle xfmr;
    status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 2.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Verify turns ratio can be retrieved */
    double ratio;
    status = mna_get_transformer_turns_ratio(&solver, xfmr, &ratio);
    ASSERT_TRUE(status == MNA_SUCCESS);
    ASSERT_DOUBLE_EQ(2.0, ratio, 1e-10);
    
    /* Change turns ratio */
    status = mna_set_transformer_turns_ratio(&solver, xfmr, 5.0);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_get_transformer_turns_ratio(&solver, xfmr, &ratio);
    ASSERT_TRUE(status == MNA_SUCCESS);
    ASSERT_DOUBLE_EQ(5.0, ratio, 1e-10);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 2: Ideal transformer validation */
static int test_ideal_transformer_validation(void) {
    MNASolver solver;
    MNAStatus status;
    ComponentHandle handle;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* Test negative turns ratio */
    status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, -2.0, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test zero turns ratio */
    status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 0.0, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test extremely large turns ratio */
    status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 1e10, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test valid turns ratio */
    status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 10.0, &handle);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    mna_destroy(&solver);
    return 1;
}

/* Run all transformer tests */
int run_transformer_tests(void) {
    printf("\n=== Transformer Tests ===\n");
    
    RUN_TEST(test_ideal_transformer_creation);
    RUN_TEST(test_ideal_transformer_validation);
    
    return tests_failed == 0;
}
