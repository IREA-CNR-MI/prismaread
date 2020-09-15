context("pr_addindex")

test_that("pr_addindex works as expected", {
    # existing index skipped ----
    testthat::expect_error(pr_addindex(Name = "MSAVI", Formula = "R600 / R700",
                                       Description = "My custom Index",
                                       Reference = "Me (2020)"))
    randind <- paste0(sprintf("%06d", sample(9999, 1, TRUE)),
                      sample(LETTERS, 1, TRUE))
    pr_addindex(Name = randind, Formula = "R600 / R700",
                Description = "My custom Index",
                Reference = "Me (2020)")

    index_list <- read.table(system.file("extdata/indexes_list.txt",
                                         package = "prismaread"), sep = "\t",
                             header = T)
    testthat::expect_true(randind %in% index_list$Name)
})

context("pr_listindexes")

test_that("pr_addindex works as expected", {
    # existing index skipped ----
   tbl <- pr_listindexes()
   testthat::expect_is(tbl, "datatables")
})

