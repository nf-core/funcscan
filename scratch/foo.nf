process printproc {
  input:
  val x

  output:
  stdout emit: hello

  script:
  """
  echo -n $x
  """
}


workflow {
  res = Channel.value([])
  res.view()
  emit:
    res
}
