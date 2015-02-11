# Test Kill function ------------------------------------------------------
n0 <- Node$new()
n1 <- n0$Spawn()
n2 <- n0$Spawn()
n21 <- n2$Spawn()
n22 <- n2$Spawn()
n0Children <- n0$GetChildren()
n2Children <- n2$GetChildren()
n2$Kill()
n0ChildrenAfterKill <- n0$GetChildren()

expect_that(length(n0ChildrenAfterKill),
            equals(length(n2Children) + length(n0Children ) - 1))