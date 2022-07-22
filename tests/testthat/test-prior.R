test_that("jsdm_prior returns right class of object", {
  expect_error(
    jsdm_prior(sigmas_b = 1),
    "All arguments must be supplied as character vectors"
  )
  expect_s3_class(jsdm_prior(), "jsdmprior")
  expect_output(print(jsdm_prior()))
})
