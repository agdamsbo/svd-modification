### Helpers

dag_shape_print <- function(dag) {
  dag |>
    ggdag::node_status() |>
    ggplot2::ggplot(ggplot2::aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = status,
      shape = status
    )) +
    ggdag:::geom_dag_edges() +
    ggdag:::geom_dag_point(size = 16) +
    ggdag:::geom_dag_label_repel(ggplot2::aes(label = label),
                                 size = 3.88, col = "black", show.legend = FALSE
    ) +
    ggplot2::labs(color="Covariate types",
                  shape="Covariate types")+
    ggdag::theme_dag()
}

coordlist2df <- function(ls){
  ls |>
    purrr::imap(\(.x, .i){
      assertthat::are_equal(length(c(.i, .x)),3)
      matrix(c(.i, .x), ncol = 3) |>
        tibble::as_tibble() |>
        setNames(c("name", "x", "y"))
    }) |> dplyr::bind_rows()
}

### DAGs

d1 <- ggdag::dagify(risk ~ pa + cvrf,
  svd ~ cvrf + risk,
  labels = c(
    "svd" = "SVD",
    "pa" = "Physical Activity",
    "cvrf" = "Cerebrovascular\nrisk-factors",
    "risk" = "CVRF risk"
  ),
  coords = list(
    "svd" = c(2, 0),
    "pa" = c(0, 2),
    "cvrf" = c(0, 0),
    risk = c(1,1)
  ) |>  coordlist2df(),
  latent = "risk",
  exposure = c("cvrf","pa"),
  outcome = "svd"
) |> ggdag::tidy_dagitty()

d2 <- ggdag::dagify(svd ~ pa,
  svd ~ cvrf,
  cvrf ~ pa,
  labels = c(
    "svd" = "SVD",
    "pa" = "Physical Activity",
    "cvrf" = "Cerebrovascular\nrisk-factors"
  ),
  coords = list(
    "svd" = c(2, 0),
    "pa" = c(1, 1),
    "cvrf" = c(0, 0)
  ) |>  coordlist2df(),
  latent = "pa",
  exposure = "cvrf",
  outcome = "svd"
) |> ggdag::tidy_dagitty()

d3 <- ggdag::dagify(svd ~ pa,
  pa ~ cvrf,
  labels = c(
    "svd" = "SVD",
    "pa" = "Physical Activity",
    "cvrf" = "Cerebrovascular\nrisk-factors"
  ),
  coords = list(
    "svd" = c(2, 0),
    "pa" = c(1, 0),
    "cvrf" = c(0, 0)
  ) |> coordlist2df(),
  latent = "pa",
  exposure = "cvrf",
  outcome = "svd"
) |> ggdag::tidy_dagitty()

### Plotting

set.seed(234)
p_dag <- list("Effect modification"=d1, 
     "Confounding"=d2, 
     "Mediation"=d3) |>
  purrr::imap(function(.x,.i){
    .x |> dag_shape_print()+
      ggplot2::labs(title=.i)+
  ggplot2::theme(legend.position = "bottom")
    }) |> 
    purrr::map2(tibble::tibble(c(-1,3),c(-1,2),c(-1,1)),function(.x,.y){
      .x+ggplot2::ylim(.y[1],.y[2])
    }) |> 
  patchwork::wrap_plots(ncol = 1,
                        guides = "collect",
                        heights = c(4,3,2,1),
                        tag_level = "keep") + 
  patchwork::guide_area()

