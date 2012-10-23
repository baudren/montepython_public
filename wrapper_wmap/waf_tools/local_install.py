def options(ctx):
  import os
  grp = ctx.parser.get_option_group("--prefix")
  ctx.parser.remove_option("--prefix")
  default_prefix = os.getcwd()
  grp.add_option("--prefix",action="store",default=default_prefix,help="installation prefix [default: %r]"%default_prefix)
  #ctx.add_option("--local",action="store_true",default=False,help="install in current directory")
  pass

def configure(ctx):
  #install where ?
  ctx.env.mprefix=ctx.env.PREFIX
  import os
  import os.path as osp
  ctx.env.localpref = os.getcwd()
  

  ctx.env.LIBDIR=osp.join(ctx.env.PREFIX,"lib")
  ctx.env.BINDIR=osp.join(ctx.env.PREFIX,"bin")
  ctx.env.INCDIR=osp.join(ctx.env.PREFIX,"include")
  
  if not os.path.exists(ctx.env.LIBDIR):
    os.mkdir(ctx.env.LIBDIR)
    
  if not os.path.exists(ctx.env.BINDIR):
    os.mkdir(ctx.env.BINDIR)
    
  if not os.path.exists(osp.join(ctx.env.PREFIX,"include")):
    os.mkdir(osp.join(ctx.env.PREFIX,"include"))
  ctx.start_msg("Setting install root to") 
  ctx.end_msg(ctx.env.PREFIX)
  ctx.start_msg("Setting install bin directory to") 
  ctx.end_msg(ctx.env.BINDIR)
  ctx.start_msg("Setting install lib directory to") 
  ctx.end_msg(ctx.env.LIBDIR)
  ctx.start_msg("Setting install include directory to") 
  ctx.end_msg(ctx.env.INCDIR)
