warn <- 0
x <- withCallingHandlers({
  devtools::document()
}, warning=function(cond) {
  warn <<-  warn+1
})

if (warn > 0) {
  stop("warnings detected in documentation")
}