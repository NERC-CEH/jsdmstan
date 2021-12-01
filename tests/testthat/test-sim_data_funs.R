test_that("gllvm_sim_data errors with bad inputs", {
  expect_error(gllvm_sim_data(N = 100, S = 10, D = 2, intercept = FALSE),
               "If K is 0 then an intercept is required")

  expect_error(gllvm_sim_data(N = "a", S = 10, D = 2),
               "N, S and D must be positive integers")
})

test_that("gllvm_sim_data returns a list of correct length", {
  expect_equal(length(gllvm_sim_data(10,2,2)), 6)
  expect_equal(length(gllvm_sim_data(10,2,2,2)), 8)
  expect_equal(length(gllvm_sim_data(10,2,2, phylo = TRUE)), 9)
  expect_equal(length(gllvm_sim_data(10,2,2,2, phylo = TRUE)), 10)
})
