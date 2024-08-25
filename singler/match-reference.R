
match_reference <- function(ref1, ref2, labels1, labels2, ...) {
  first <- singler(test = ref1, ref = ref2, labels = labels2, ...)
  second <- singler(test = ref2, ref = ref1, labels = labels1, ...)

  f1 <- factor(labels1)
  f2 <- factor(labels2)
  tab1 <- table(f1, factor(first$labels, levels(f2)))
  tab2 <- table(f2, factor(second$labels, levels(f1)))

  tab1 <- tab1 / rowSums(tab1)
  tab2 <- tab2 / rowSums(tab2)
  output <- tab1 * t(tab2)

  output <- unclass(output)
  dimnames(output) <- unname(dimnames(output))
  output
}
