dag <- dagitty::dagitty('dag {
  A [selected,pos="-2.200,-1.520"]
  B [pos="1.400,-1.460"]
  D [outcome,pos="1.400,1.621"]
  E [exposure,pos="-2.200,1.597"]
  Z [adjusted,pos="-0.300,-0.082"]
  A -> E
  A -> Z [pos="-0.791,-1.045"]
  B -> D
  B -> Z [pos="0.680,-0.496"]
  E -> D
}')

tidy_dag <- ggdag::tidy_dagitty(dag)

ggdag::ggdag(tidy_dag) +
  ggdag::theme_dag()

