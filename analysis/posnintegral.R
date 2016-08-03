require(dplyr)
# Do "Integral plots" of points to show range
position_integral <- function(width = NULL, height = NULL) {
  ggproto(NULL, PositionIntegral,
          width = width,
          height = height
  )
}

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

PositionIntegral <- ggproto("PositionIntegral", Position,
                          required_aes = c("x", "y"),
                          
                          setup_params = function(self, data) {
                            list(
                              width = self$width %||% resolution(data$x, zero = FALSE) * 0.4,
                              height = self$height %||% resolution(data$y, zero = FALSE) * 0.4
                            )
                          },
                          
                          compute_layer = function(data, params, panel) {
                            # I don't really understand why the arrange(y) is necessary here
                            data <- data %>% 
                                      group_by(x,PANEL) %>% 
                                      arrange(y) %>%
                                      mutate(o=order(y), dx=(o-1)/(max(o)-1)) %>% 
                                      mutate(newx = x + (dx-0.5)*2*params$width ) %>%
                                      ungroup()
                            data$x <- data$newx
                            data
                          }
)