test_that("Joining tables from nt_descriptive by group and overall", {
  withr::local_package("dplyr")

  iris_nt <- iris |> filter(Species != "versicolor") |>
    mutate(Species = droplevels(Species))
  tab01 <- nt_describe(iris_nt, group = Species)
  tab02 <- nt_describe(iris_nt)
  tab <- nt_join_tables(tab01, tab02)
  expect_snapshot(as.data.frame(tab))
})

test_that("Joining tables from nt_descriptive and nt_compare_tg", {
  withr::local_package("dplyr")

  iris_nt <- iris |> filter(Species != "versicolor") |>
    mutate(Species = droplevels(Species))
  tab01 <- nt_describe(iris_nt, group = Species)
  tab02 <- nt_compare_tg(iris_nt, group = Species)
  tab <- nt_join_tables(tab01, tab02)
  expect_snapshot(as.data.frame(tab))
})

test_that("Joining tables from nt_descriptive and nt_compare_mg", {
  withr::local_package("dplyr")

  iris_nt <- iris
  tab01 <- nt_describe(iris_nt, group = Species)
  tab02 <- nt_compare_mg(iris_nt, group = Species)
  tab <- nt_join_tables(tab01, tab02)
  expect_snapshot(as.data.frame(tab))
})
