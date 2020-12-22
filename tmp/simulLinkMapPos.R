
simulLinkMapPos <- function(pos, len, lenCM = 100, b1 = 10, b2 = 10){
  stopifnot(b1 >= 2)
  stopifnot(b2 >= 2)

  posCM <- ((pos/len + b1 #- len/2
             - dbeta((pos/len), 1, b1)
             + dbeta((pos/len), b2, 1))
            * lenCM/(1 + b1 + b2))
  return(posCM)
}



len <- 1000
lenCM <- 100
n <- 10
b1 = 10
b2 = 10000
SNPcoord <- data.frame(physPos = sort(sample(len, size = n, replace = FALSE)))
SNPcoord$linkMapPos <- simulLinkMapPos(SNPcoord$physPos,
                                       len,
                                       lenCM,
                                       b1 = b1,
                                       b2 = b2)

is.unsorted(SNPcoord$linkMapPos, na.rm = FALSE, strictly = FALSE)



(nRecomb <- rpois(1, lenCM/100))

RposCm <- runif(nRecomb, 0, lenCM)
# get id of the recombination position
ids <- (c(0,
          sort(base::findInterval(RposCm, SNPcoord$linkMapPos)),
          nrow(SNPcoord)))

SNPcoord$grp <- 0

g <- 1
for (i in seq_len(length(ids) - 1)) {
  if (ids[i] == ids[i + 1]) {
    g <- -g + 3
    next
  }
  SNPcoord$grp[(ids[i] + 1):ids[i + 1]] <- g
  g <- -g + 3
}



plotDta <- data.frame(x = seq(0, len, length.out = 200)) %>%
  dplyr::mutate(y1 = simulLinkMapPos(x, len,lenCM,b1 = b1,b2 = b2))
p <- plot_ly(
  type = "scatter",
  mode = "lines",
  data = plotDta,
  x = ~ x,
  y = ~ y1,
  name = "y1"
)  %>%
  add_markers(data = SNPcoord,
              x = ~physPos,
              y= ~linkMapPos,
              color = ~as.character(grp)) %>%
  # add_segments(x = 0, xend = len, y = RposCm, yend = RposCm) %>%
  layout(hovermode = 'compare')
p
