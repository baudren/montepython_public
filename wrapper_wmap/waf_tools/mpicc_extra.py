def configure(ctx):
  ctx.env.has_mpi = False
  mpiccpath = ctx.find_program("mpicc")
  if mpiccpath:
    ctx.env.has_mpi = True
    envmpi = ctx.env.copy() 
    ctx.setenv('mpi', envmpi) 
    ctx.env.CC = [mpiccpath]
    ctx.env.LINK_CC = [mpiccpath]
    envmpibld = envmpi = ctx.env.copy() 
    ctx.set_env_name('mpibld', envmpibld) 
