source(here::here("R/functions.R"))

targets::tar_read("df_complete") |> 
  ggplot2::ggplot(ggplot2::aes(x=month,y=pase_0))+
  ggplot2::geom_boxplot()+
  ggplot2::geom_point()

targets::tar_read("df_complete") |> 
  ggplot2::ggplot(ggplot2::aes(x=winter,y=pase_0))+
  ggplot2::geom_boxplot()+
  ggplot2::geom_point()


targets::tar_read("df_complete") |> 
  (\(.x){
    t.test(pase_0~winter,data=.x)
  })()
  
  


targets::tar_read("df_complete") |> 
  ggplot2::ggplot(ggplot2::aes(x=year,
                               y=pase_0, 
                               # color=trial, 
                               fill=active_treatment))+
  ggplot2::geom_boxplot()+
  ggplot2::geom_point()
