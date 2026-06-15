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

df1 <- DataFrame(a=x1, b=x2, c=x3)
df2 <- df1
df2$d <- y1
df2$e <- y2
df2$f <- y3

test_that(".select_consistent_columns works correctly", {
  expect_equal(.select_consistent_columns(df2, rv), df1)
})

df <- data.frame(a=c('a', 'b', 'c', 'd'),
                 b=c('a', 'a', 'a', 'b'),
                 c=c('d', 'd', 'd', 'd'),
                 d=c('a', 'a', 'a', 'b'),
                 e=c('b', 'b', 'a', 'a'))

test_that(".select_consistent_columns works correctly for vectors whose factorization has a lower resolution than the cluster factorization",
          {
              expect_identical(.select_consistent_columns(df, df$a), df)
              expect_identical(.select_consistent_columns(df, df$e), df[,c("c", "e")])
          })



df1 <- DataFrame(a=y1, b=x2, c=y3)
df1$d <- DataFrame(g=x2, f=x3, h=y2)
df2 <- DataFrame(b=x2)
df2$d <- DataFrame(g=x2, f=x3)

test_that(".select_consistent_columns works correctly with nested DataFrame", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})


df1 <- DataFrame(a=y1, b=x2, c=y3)
df1$d <- DataFrame(g=x2)
df1$e <- DataFrame(g=y2)
df2 <- DataFrame(b=x2)
df2$d <- DataFrame(g=x2)

test_that(".select_consistent_columns works correctly with nested DataFrames composed of the one column", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})

df1 <- DataFrame(a=y1, b=x2, c=y3)
df1$d <- matrix(c(y2,x1), ncol=2)
df2 <- DataFrame(b=x2)
df2$d <- matrix(x1, ncol=1)

df3 <- S4Vectors::DataFrame('a'=c('a','a','a','b'))
df3$b <- matrix(rep(c(c('c','c','c','a')),3), nrow =4)
df3[["c"]] <- cbind(df3$b, c('a', 'b', 'c', 'd'))
df4 <- df3
df4[[3]] <- df4[[3]][,1:3]

test_that(".select_consistent_columns works correctly with nested matrix", {
    expect_equal(.select_consistent_columns(df1, rv), df2)
    expect_identical(.select_consistent_columns(df3, df3$a), df4)
    expect_identical(.select_consistent_columns(df3, c('a', 'b', 'c', 'd')), df3)
    expect_identical(.select_consistent_columns(df3,  rep('a',4)), NULL)
})

df1 <- DataFrame(a=y1, b=x2, c=y3)
df1$d <- matrix(x2, ncol=1)
df1$e <- matrix(y2, ncol=1)
df2 <- DataFrame(b=x2)
df2$d <- matrix(x2, ncol=1)

test_that(".select_consistent_columns works correctly with nested matrix composed of the one column", {
  expect_equal(.select_consistent_columns(df1, rv), df2)
})
