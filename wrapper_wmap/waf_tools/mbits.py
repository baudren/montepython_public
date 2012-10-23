def options(ctx):

  import optparse
  grp=optparse.OptionGroup(ctx.parser,"architecture options")
  grp.add_option("--m32",action="store_true",default=False,help="compile & link in 32bits")
  grp.add_option("--m64",action="store_true",default=False,help="compile & link in 64bits [default]")
  ctx.add_option_group(grp)  


def configure(ctx):
  import sys
  from waflib import Utils,Errors
  if ctx.options.m32 and ctx.options.m64:
    raise Errors.WafError("You must choose either m32 of m64 !")
  #32 bits
  if ctx.options.m32==False and ctx.options.m64==False:
    ctx.options.m64=True
  if sys.platform.lower()=="darwin":
    mopt = ""
    if ctx.options.m64:
      mopt += "-arch x86_64 "
    if ctx.options.m32:    
      mopt += "-arch i386 "    
  else:
    mopt = "-m64"
    if ctx.options.m32:
      mopt = "-m32"
      
  ctx.env.mopt=mopt
  ctx.env.append_value('CCFLAGS',mopt)
  ctx.env.append_value('LINKFLAGS',mopt)
  ctx.start_msg("Setting architecture flag to") 
  ctx.end_msg(mopt)
