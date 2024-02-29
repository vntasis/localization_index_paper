#' npidr
#'
#' npIDR [non-parametric Irreproducible Discovery Rate]
#'

npidr <- function(input, bin.size = 1) {
  data <- round(input / bin.size, digits = 0)

  acount1 <- as.data.frame(table(data[, 1]))
  ccount1 <- as.data.frame(table(data[data[, 2] == 0, 1]))

  matr1 <- merge(acount1, ccount1, by = 1, all = TRUE)

  matr1[is.na(matr1)] <- 0

  matr1$idr <- matr1$Freq.y / matr1$Freq.x

  acount2 <- as.data.frame(table(data[, 2]))
  ccount2 <- as.data.frame(table(data[data[, 1] == 0, 2]))

  matr2 <- merge(acount2, ccount2, by = 1, all = TRUE)

  matr2[is.na(matr2)] <- 0

  matr2$idr <- matr2$Freq.y / matr2$Freq.x

  npidr1 <- matr1$idr
  names(npidr1) <- matr1[, 1]

  npidr2 <- matr2$idr
  names(npidr2) <- matr2[, 1]

  npidr <- (npidr1[as.character(data[, 1])]
            + npidr2[as.character(data[, 2])]) / 2
  npidr <- as.numeric(npidr)
  return(npidr)
}
