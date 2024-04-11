# reference vector
rv <- c("a", "b", "b", "a", "c")

# consistent vectors
x1 <- rv
x2 <- c("b", "a", "a", "b", "c")
x3 <- c("b", "c", "c", "b", "c")

# inconsistent vectors
y1 <- c("v", "b", "b", "a", "c")
y2 <- c("b", "b", "b", "a", "c")
y3 <- c("a", "b", "a", "a", "a")

df1 <- DataFrame(a = x1, b = x2, c = x3)
df2 <- df1
df2$d <- y1
df2$e <- y2
df2$f <- y3

test_that(".select_consistent_columns works correctly", {
  expect_equal(.select_consistent_columns(df2, rv), df1)
})

df1 <- DataFrame(a = y1, b = x2, c = y3)
df1$d <- DataFrame(g = x2, f = x3, h = y2)
df2 <- DataFrame(b = x2)
df2$d <- DataFrame(g = x2, f = x3)

test_that(".select_consistent_columns works correctly with nested DataFrame", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})


df1 <- DataFrame(a = y1, b = x2, c = y3)
df1$d <- DataFrame(g = x2)
df1$e <- DataFrame(g = y2)
df2 <- DataFrame(b = x2)
df2$d <- DataFrame(g = x2)

test_that(".select_consistent_columns works correctly with nested DataFrames composed of the one column", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})

df1 <- DataFrame(a = y1, b = x2, c = y3)
df1$d <- matrix(c(y2,x1), ncol=2)
df2 <- DataFrame(b = x2)
df2$d <- matrix(x1, ncol=1)

test_that(".select_consistent_columns works correctly with nested matrix", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})

df1 <- DataFrame(a = y1, b = x2, c = y3)
df1$d <- matrix(x2, ncol=1)
df1$e <- matrix(y2, ncol=1)
df2 <- DataFrame(b = x2)
df2$d <- matrix(x2, ncol=1)

test_that(".select_consistent_columns works correctly with nested matrix composed of the one column", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})
