import re
import sys
from waflib import Options
import os.path as osp
from waflib import Logs
from waflib import  Context
from waflib import Errors

def options(ctx):
  import optparse
  grp = ctx.parser.get_option_group("--gcc")
  if grp==None:
    grp=optparse.OptionGroup(ctx.parser,"compiler options")
  grp.add_option("--gfortran",action="store_true",default=False,help="Do not test for ifort and only use gfortran")
  grp.add_option("--ifort",action="store_true",default=False,help="Do not test for gfortran and only use ifort")
  grp.add_option("--fortran_flagline",action="store",default="",help="flagline to link fortran object to c using ld")
  ctx.add_option_group(grp)  
  
def configure_(ctx):
  if ctx.options.fortran_flagline:
    conf.parse_flags(ctx.options.fortran_flagline,uselib="fc_runtime")
  if sys.platform.lower()=="darwin":
    ctx.env.fcshlib_PATTERN = 'lib%s.dylib'
  
  
  ctx.env.has_ifort = False
  if not Options.options.gfortran:
    try:
      ifort_conf(ctx)
      return
    except Exception,e:
      if Options.options.ifort:
        raise
      Logs.pprint("PINK", "ifort not found, defaulting to gfortran (cause: '%s')"%e)
  gfortran_conf(ctx)

def configure(ctx): 
  configure_(ctx) 
  ctx.env.append_value("FCFLAGS_fcshlib",ctx.env.LINKFLAGS_fcshlib)  
  ctx.env["FCFLAGS_fpic"]=[]
  ctx.env.append_value("FCFLAGS_fpic",[flg for flg in ctx.env.FCFLAGS_fcshlib if "-fpic" in flg.lower()])
  
def show_linkline(ctx):
  ctx.start_msg("fortran link line")
  ctx.end_msg(" ".join(["-L%s"%vv for vv in ctx.env.LIBPATH_fc_runtime])+" "+" ".join(["-l%s"%vv for vv in ctx.env.LIB_fc_runtime]))
  
def ifort_conf(ctx):
  import waflib
  import os
  ctx.env.FC=[]
  ctx.load('ifort')
  if sys.platform.lower()=="darwin":
    ctx.env.LINKFLAGS_fcshlib = ['-dynamiclib']
  ctx.env.append_value('FCFLAGS',ctx.env.mopt.split())
  ctx.env["FCFLAGS_fc_omp"]=[]
  ctx.env.append_value("FCFLAGS_fc_omp","-openmp")
  ctx.env.FCSHLIB_MARKER = [""]
  ctx.env.FCSTLIB_MARKER = [""]
  ctx.check_cc(
    errmsg="failed",msg='Compile a test code with ifort',
    mandatory=1,fragment = "program test\n  WRITE(*,*) 'hello world'\n end program test\n",compile_filename='test.f90',features='fc fcprogram')
  if not ctx.options.fortran_flagline:
    ctx.start_msg("retrieve ifort link line")
    try:
      #print "%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.f90"%(ctx.env.FC," ".join(ctx.env.FCFLAGS))
      llgo,llge = ctx.cmd_and_log("%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.f90"%(ctx.env.FC," ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fc_omp)), output=waflib.Context.BOTH)
      #print "RET",llgo,llge
      L = set([ll.strip() for ll in re.findall("-L(.+)\s*\\\\", llge.split("ld ")[1]) if ("ifort" in ll.lower()) or ("intel" in ll.lower())])
      l = set([ll.strip() for ll in re.findall("-l(.+)\s*\\\\", llge.split("ld ")[1])])
      rL = set()
      rl = set()
      for Li in L:
        oli = os.listdir(Li)
        for li in l:
          if ctx.env.cshlib_PATTERN%li in oli:
            rl.add(li)
            rL.add(Li)
    except:
      ctx.end_msg(False)
      raise
    for pth in list(rL) + ["/lib","/lib64"]:
      if osp.exists(pth):
        ctx.env.append_value("LIBPATH_fc_runtime",pth)
        ctx.env.append_value("RPATH_fc_runtime",pth)
    
    ctx.env.append_value("LIB_fc_runtime",list(rl)+["pthread"])
    ctx.end_msg(True)
  show_linkline(ctx)
  ctx.env.has_ifort = True

def ifort_conf_(ctx):
  ctx.env.FC=[]
  ctx.load('ifort')
  if sys.platform.lower()=="darwin":
    ctx.env.LINKFLAGS_fcshlib = ['-dynamiclib']
  ctx.env.append_value('FCFLAGS',ctx.env.mopt.split())
  ctx.env.append_value("FCFLAGS_fc_omp","-openmp")
  ctx.env.FCSHLIB_MARKER = [""]
  ctx.env.FCSTLIB_MARKER = [""]
  ctx.check_cc(
    errmsg="failed",msg='Compile a test code with ifort',
    mandatory=1,fragment = "program test\n  WRITE(*,*) 'hello world'\n end program test\n",compile_filename='test.f90',features='fc fcprogram')
  if not ctx.options.fortran_flagline:
    ctx.start_msg("retrieve ifort link line")
    if "/" not in ctx.env.FC:
      ctx.env.FC = ctx.cmd_and_log("which %s"%ctx.env.FC).strip()
      #print ctx.env.FC
    ifort_path = osp.dirname(osp.realpath(ctx.env.FC))
    
    #print ifort_path
    if ctx.options.m32:
      try:
        f=open(osp.join(ifort_path,'ifortvars_ia32.sh'))
      except:
        ctx.end_msg(False)
        raise Errors.WafError("Can't locate ifort configuration file")
    else:
      try:
        f=open(osp.join(ifort_path,'ifortvars_intel64.sh'))
      except:
        ctx.end_msg(False)
        raise Errors.WafError("Can't locate ifort configuration file")

    txt = f.read()
    f.close()
    #print txt
    if sys.platform.lower()=="darwin":
      sp = "DYLD_LIBRARY_PATH"
    else:
      sp = "LD_LIBRARY_PATH"
    res = re.findall("\s"+sp+"\s*=\s*\"(.+)\"",txt)[0]
    for pth in res.split(":"):
      ctx.env.append_value("LIBPATH_fc_runtime",pth)
      ctx.env.append_value("RPATH_fc_runtime",pth)
    ctx.env.append_value("LIB_fc_runtime",["ifcore","intlc","ifport","imf","irc","svml","iomp5","pthread"])
    ctx.end_msg(True)
  show_linkline(ctx)  
  
def gfortran_conf(ctx):
  ctx.env.FC=[]
  ctx.env.FCFLAGS = []
  ctx.load('gfortran')
  ctx.env["FCFLAGS_fc_omp"]=[]
  ctx.env.append_value("FCFLAGS_fc_omp","-fopenmp")
  ctx.env.append_value("FCFLAGS","-DGFORTRAN")
  ctx.env.append_value("FCFLAGS","-ffixed-line-length-0")
  ctx.env.append_value("FCFLAGS","-ffree-line-length-0")
  mopt = ctx.env.mopt
  if sys.platform.lower()=="darwin":
    if "i386" in ctx.env.mopt:
      ctx.env.append_value('FCFLAGS','-m32')
      mopt = "-m32"
    else:
      ctx.env.append_value('FCFLAGS','-m64')
      mopt = "-m64"
  else:
    ctx.env.append_value('FCFLAGS',ctx.env.mopt.split())
  ctx.start_msg("Check gfortran version") 
  v90 = ctx.cmd_and_log(ctx.env.FC+" --version",quiet=Context.STDOUT).split("\n")[0].strip()
  version90 = re.findall("(4\.[0-9]\.[0-9])",v90)
  if len(version90)<1:
    #Logs.pprint("PINK","Can't get gfortran version... Let's hope for the best")
    ctx.end_msg("not found, let's hope for the best...",color="PINK")
  else:
    version90 = version90[0]
    vmid = int(version90.split(".")[1])
    if vmid<3:
      ctx.end_msg(v90,color="YELLOW")
      raise Errors.WafError("gfortran version need to be above 4.3 got %s"%version90)
    ctx.end_msg(v90)
  
  # kludge !
  ctx.env.FCSHLIB_MARKER = [""]
  ctx.env.FCSTLIB_MARKER = [mopt]
  ctx.check_cc(
      errmsg="failed",msg='Compile a test code with gfortran',
      mandatory=1,fragment = "program test\n  WRITE(*,*) 'hello world'\n end program test\n",compile_filename='test.f90',features='fc fcprogram')

  ctx.start_msg("retrieve gfortran link line")
  lgfpath = ctx.cmd_and_log(ctx.env.FC+" %s -print-file-name=libgfortran.dylib"%mopt,quiet=Context.STDOUT)    
  lpath = [osp.dirname(osp.realpath(lgfpath))]
  lgfpath = ctx.cmd_and_log(ctx.env.FC+" %s -print-file-name=libgomp.dylib"%mopt,quiet=Context.STDOUT)    
  lpath += [osp.dirname(osp.realpath(lgfpath))]
  lpath = set(lpath)

  ctx.env.append_value("LIB_fc_runtime",["gfortran","gomp"])
  ctx.env.append_value("LIBPATH_fc_runtime",list(lpath))
  ctx.env.append_value("RPATH_fc_runtime",list(lpath))
  ctx.end_msg(True)
  
  show_linkline(ctx)