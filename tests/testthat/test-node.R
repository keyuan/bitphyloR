NodeTestTemplete <- function(f) {
  test_that("Node functions tests", {
    test_that("classname contains 'Node'", {
      n0 <- f$new()
      expect_that("Node" %in% class(n0), is_true())
    })

    test_that("AddChild function adds a child node",{
      n0 <- f$new()
      n1 <- f$new()
      n0$AddChild(n1)
      expect_that(n0$GetChildren()[[1]], is_identical_to(n1))
      expect_that(n1$GetParent(), is_identical_to(n0))
    })

    test_that("Remove a child node", {
      n0 <- f$new()
      n1 <- f$new(n0)
      expect_that(n0$GetChildren()[[1]], is_identical_to(n1))
      n0$RemoveChild(n1)
      expect_that(length(n0$GetChildren()), is_equivalent_to(0))
    })

    test_that("Spawn function generates a Node child", {
      n0 <- f$new()
      n1 <- n0$Spawn()
      n0Children <- n0$GetChildren()
      expect_that("Node" %in% class(n1), is_true())
      expect_that(length(n0Children), equals(1))
      expect_that(n0$GetChildren()[[1]], is_identical_to(n1))
    })

    test_that("Kill function kills and pushes children",{
      n0 <- f$new()
      n1 <- n0$Spawn()
      n1$Kill()
      expect_that(length(n0$GetChildren()), is_equivalent_to(0))
    })
  })
}

testObjects <- c(Node, Normal)
lapply(testObjects, NodeTestTemplete)
