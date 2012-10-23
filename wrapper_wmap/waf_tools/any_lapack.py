#try to support many flavours of lapack
import autoinstall_lib as atl
from waflib import Logs
import os.path as osp
from waflib import Logs
from waflib import  Context
from waflib import Errors
import sys
import waflib

version = "lapack-3.3.1"
tool = "lapack-3.3.1"
lapack_funcs = "dtrsv dpotrf dpotrs dpotri dtrtri dtrmm dtrmv dgeqrf dormqr dsyev dgesvd dsymv dgemv dgemm dsyrk dsyr2k daxpy dtrsm dsymm dsyr ddot"     

def options(ctx):
  atl.add_lib_option("lapack",ctx,install=True)
  grp = ctx.parser.get_option_group("--lapack_install")
  grp.add_option("--lapack_mkl",action="store",default="",help="if lapack is mkl, location of the mkl install")
  grp.add_option("--lapack_mkl_version",action="store",default="",help="version of the mkl lib (should be 10.3, 10.2, 10.1 or 10.0)")
  grp.add_option("--lapack_apple",action="store_true",default=False,help="use apple version of blas/lapack")

def do_include(ctx,ptrn="%s_"):
  f=open(osp.join(ctx.env.PREFIX,"include/lapack_clik.h"),"w")
  for fnc in lapack_funcs.split():
    print >>f,("#define %s "+ptrn)%(fnc,fnc)
  print >>f,extra_inc
  f.close()

def configure(ctx):
  #always assume that I need a dedicated include file.
  
  if ctx.options.lapack_apple:
    ctx.start_msg("Check apple lapack")
    if sys.platform.lower()!="darwin":
      ctx.end_msg("not on darwin ! Got '%s'"%sys.platform,color="YELLOW")
      raise Errors.WafError("cannot find apple lapack")
    ctx.end_msg("ok")
    lapack_extradefs = ["HAS_LAPACK"]
    lapack_libs = ["BLAS","LAPACK"]
    lapack_includes = ["lapack_clik.h"]
    lapack_extradefs += ["LAPACK_CLIK"]
    ctx.options.lapack_include = osp.join(ctx.env.PREFIX,"include")
    ctx.options.lapack_lib = "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current"
    do_include(ctx,"%s_")

  elif ctx.options.lapack_mkl:
    # parse version
    ctx.start_msg("Check mkl version")
    if ctx.options.lapack_mkl_version.strip()[:4] not in ("10.0","10.1","10.2","10.3"):
      ctx.end_msg(ctx.options.lapack_mkl_version.strip(),color="YELLOW")
      raise Errors.WafError("Cannot understand mkl version '%s'"%ctx.options.lapack_mkl_version.strip())
    version = int(ctx.options.lapack_mkl_version.strip()[:4].split(".")[1])
    ctx.end_msg("10.%d"%version)
    lapack_extradefs = ["HAS_LAPACK"]
    lapack_extradefs += ["HAS_MKL"]
    lapack_includes = ["mkl_lapack.h","mkl_blas.h"]
    lapack_libs = []
    tag = sys.platform.lower()
    if tag=="darwin":
      pass
    elif "linux" in tag:
      tag="linux"
    else:
      raise Errors.WafError("unknown platform '%s'"%tag)
    tag+="_10.%d"%version
    mopt = ctx.env.mopt
    if "64" in mopt:
      tag+="_64"
    else:
      tag +="_32"

    if sys.platform.lower()!='darwin':
      #I need to create my own lapack ! 
      cmdline = """gcc -shared -Bdynamic  %(func_list)s  -Wl,--start-group %(ars)s  -Wl,--end-group %(Lomp)s %(omp)s -o "%(res)s" """
      cmdlist = {}
      cmdlist["func_list"] = " ".join(["-u %s_"%v for v in lapack_funcs.split()])
      cmdlist["ars"] = " ".join([osp.join(mkl_options[tag][0]%(ctx.options.lapack_mkl),"lib%s.a"%v.strip()) for v in mkl_options[tag][1].split("-l") if v.strip() and v.strip()[:3]=="mkl"])
      cmdlist["Lomp"] = " ".join("-L%s"%v.strip() for v in ctx.env.LIBPATH_fc_runtime if v.strip())
      cmdlist["omp"] = " ".join([v.strip() for v in mkl_options[tag][1].split() if v.strip() and "mkl" not in v])
      cmdlist["res"] = osp.join(ctx.env.LIBDIR,ctx.env.cshlib_PATTERN%"clik_mkl")
      cmdline = cmdline%cmdlist
      #print cmdline
      ctx.start_msg("create specific mkl lib")
      llgo,llge = ctx.cmd_and_log(cmdline, output=waflib.Context.BOTH)
      #print llgo
      #print llge
      ctx.end_msg(cmdlist["res"])
      ctx.options.lapack_link = "-lclik_mkl "+cmdlist["omp"]
      ctx.options.lapack_lib = ctx.env.LIBDIR+":".join([""]+ctx.env.LIBPATH_fc_runtime)
      ctx.options.lapack_include =  ctx.options.lapack_mkl+"/include"

    else:
      ctx.options.lapack_link = mkl_options[tag][1]
      ctx.options.lapack_lib = mkl_options[tag][0]%(ctx.options.lapack_mkl)+":".join([""]+ctx.env.LIBPATH_fc_runtime)
      if "framework" in ctx.options.lapack_mkl.lower():
        ctx.options.lapack_include =  ctx.options.lapack_mkl+"/Headers"
      else:
        ctx.options.lapack_include =  ctx.options.lapack_mkl+"/include"

    #try:
    #  atl.conf_lib(ctx,"lapack",lapack_libs,lapack_funcs.split(),lapack_includes,defines=lapack_extradefs,install=installlapack)
    #except Exception,e:
    #  pass

  #lapack_extradefs = ["HAS_LAPACK"]
  #lapack_libs = ["BLAS","LAPACK"]
  #lapack_includes = ["lapack.h","blas.h"]

  #if "mkl" in ctx.options.lapack_lib.lower() or "mkl" in ctx.options.lapack_include.lower() or "mkl" in ctx.options.lapack_link or ctx.options.lapack_mkl:
  #  ctx.env.mkl = True
  #  lapack_extradefs += ["HAS_MKL"]
  #  lapack_includes = ["mkl_lapack.h","mkl_blas.h"]
  #  if ctx.options.lapack_mkl:
  #    if ctx.env.has_ifort==False:
  #      raise Exception("cannot use MKL without ifort")
  #    if "framework" in ctx.options.lapack_mkl.lower():
  #      # guess we are on macosx
  #      # get the path of the framework
  #      if ctx.options.lapack_mkl[-1] == "/":
  #        fpath,fname = osp.split(ctx.options.lapack_mkl[:-1])
  #      else:
  #        fpath,fname = osp.split(ctx.options.lapack_mkl)
  #      fname = fname.split(".")[0]
  #      ctx.options.lapack_include =  ctx.options.lapack_mkl+"/Headers"
  #      ctx.options.lapack_lib =  ctx.options.lapack_mkl+"/Libraries/universal"
  #      if ctx.options.lapack_link=="":
  #        ctx.options.lapack_link = "-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core"
  #    else:
  #      # assume it's 10 on linux
  #      # check whether it's 10.3
  #      if ctx.options.m32:
  #        libsuffix="/lib/32"
  #        libdep = "-lmkl_intel"
  #      else:
  #        libsuffix="/lib/em64t"
  #        libdep = "-lmkl_intel_lp64"
  #      if ctx.options.lapack_link=="":
  #        ctx.options.lapack_link = "-lmkl_lapack -lmkl_intel_thread -lmkl_core -liomp5 -lm -lpthread -lmkl_def" + libdep
  #      if not ctx.options.m32 and osp.exists(ctx.options.lapack_mkl+"/lib/intel64"):
  #        libsuffix="/lib/intel64"
  #        ctx.options.lapack_link = "-lmkl_intel_thread -lmkl_core -liomp5 -lm -lpthread -lmkl_def" + libdep
  #      ctx.options.lapack_include=ctx.options.lapack_mkl+"/include"
  #      ctx.options.lapack_lib=ctx.options.lapack_mkl+libsuffix+":".join([""]+ctx.env.LIBPATH_fc_runtime)
  elif atl.upgrade(ctx,"lapack") or ctx.options.lapack_islocal or ctx.options.lapack_forceinstall or atl.shouldIinstall_all(ctx,"lapack"):
    ctx.env.append_value("LIBPATH_lapack",ctx.env.LIBPATH_fc_runtime)
    ctx.env.append_value("RPATH_lapack",ctx.env.RPATH_fc_runtime)
    ctx.env.append_value("LIB_lapack",ctx.env.LIB_fc_runtime)
    lapack_libs = ["lapack_clik","blas_clik"]
    lapack_includes = ["lapack_clik.h"]
    lapack_extradefs = ["HAS_LAPACK"]
    lapack_extradefs += ["LAPACK_CLIK"]
  else:
    lapack_libs = []
    lapack_includes = ["lapack_clik.h"]
    lapack_extradefs = ["HAS_LAPACK"]
    lapack_extradefs += ["LAPACK_CLIK"]
    do_include(ctx)

  atl.conf_lib(ctx,"lapack",lapack_libs,lapack_funcs.split(),lapack_includes,defines=lapack_extradefs,install=installlapack)

def installlapack(ctx):
  filen = version+".tgz"
  atl.installsmthg_pre(ctx,"http://www.netlib.org/lapack/"+filen,filen)
  from waflib import Utils,Errors
  dii = {"FCC":ctx.env.FC,"FCFLAGS":" ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fcshlib),"FLINKFLAGS":" ".join(ctx.env.FCFLAGS+ctx.env.LINKFLAGS_fcshlib),"SO":ctx.env.shsuffix,"MFLAG":" ".join(ctx.env.FCFLAGS) }
  Logs.pprint("PINK","build blas")
  f=open("build/%s/make.inc"%version,"w")
  print >>f,make_inc_blas%dii
  f.close()
  cmdline = "cd build/%s; make blaslib"%version
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%version)
  Logs.pprint("PINK","build lapack")
  f=open("build/%s/make.inc"%version,"w")
  print >>f,make_inc_lapack%dii
  f.close()
  cmdline = "cd build/%s; make lapacklib"%version
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%version)
  
  import shutil
  shutil.copyfile("build/%s/liblapack_clik.%s"%(version,ctx.env.shsuffix), osp.join(ctx.env.LIBDIR,"liblapack_clik.%s"%ctx.env.shsuffix))
  shutil.copyfile("build/%s/libblas_clik.%s"%(version,ctx.env.shsuffix), osp.join(ctx.env.LIBDIR,"libblas_clik.%s"%ctx.env.shsuffix))

  do_include(ctx)
  
make_inc_lapack="""
SHELL = /bin/sh
FORTRAN  = %(FCC)s %(FCFLAGS)s
OPTS     =
DRVOPTS  = $(OPTS)
NOOPT    = -g -O0
TIMER    = INT_CPU_TIME
LOADER   = %(FCC)s
LOADOPTS = %(MFLAG)s

BLASLIB      = ../../libblas_clik.%(SO)s
ARCH = %(FCC)s 
ARCHFLAGS = %(FLINKFLAGS)s -L../ -lblas_clik -o
RANLIB = echo
LAPACKLIB    = liblapack_clik.%(SO)s
"""

make_inc_blas="""
SHELL = /bin/sh
FORTRAN  = %(FCC)s %(FCFLAGS)s
OPTS     =
DRVOPTS  = $(OPTS)
NOOPT    = -g -O0
TIMER    = INT_CPU_TIME

BLASLIB      = ../../libblas_clik.%(SO)s
ARCH = %(FCC)s 
ARCHFLAGS = %(FLINKFLAGS)s -o
RANLIB = echo
LAPACKLIB    = liblapack_clik.%(SO)s
"""

extra_inc = """
void dtrsv(const char *uplo, const char *trans, const char *diag, const int  *n,
           const double *a, const int *lda, double *x, const int *incx);
void dpotrf( char* uplo, int * n, double* a, int * lda, int * info );
void dpotri( char* uplo, int * n, double* a, int * lda, int * info );
void dgemv(const char *trans, const int *m, const int *n, const double *alpha,
           const double *a, const int *lda, const double *x, const int *incx,
           const double *beta, double *y, const int *incy);
void dsyrk(const char *uplo, const char *trans, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *beta,
           double *c, const int *ldc);
void dsyr2k(const char *uplo, const char *trans, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);
void dgesvd( char* jobu, char* jobvt, int * m, int * n, double* a, int * lda, double* s, double* u, int * ldu, double* vt, int * ldvt, double* work, int * lwork, int * info );
void dgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);
void dtrtri( char* uplo, char* diag, int * n, double* a, int * lda, int * info );
void dtrmm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);
void dtrmv(const char *uplo, const char *transa, const char *diag, const int *n,
           const double *a, const int *lda, double *b, const int *incx);
void dgeqrf( int * m, int * n, double* a, int * lda, double* tau, double* work, int * lwork, int * info );
void dormqr( char* side, char* trans, int * m, int * n, int * k, double* a, int * lda, double* tau, double* c, int * ldc, double* work, int * lwork, int * info );
void dsyev( char* jobz, char* uplo, int * n, double* a, int * lda, double* w, double* work, int * lwork, int * info );
void dsymv(const char *uplo, const int *n, const double *alpha, const double *a, const int *lda,
           const double *x, const int *incx, const double *beta, double *y, const int *incy);
void daxpy(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
void dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const double *alpha, const double *a, const int *lda,
           double *b, const int *ldb);
void dsyr(const char *uplo, const int *n, const double *alpha, const double *x, const int *incx,
         double *a, const int *lda);
void dsymm(const char *side, const char *uplo, const int *m, const int *n,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);
double ddot(int* N,double *DX, int* INCX,double *DY,int* INCY);
void dpotrs(char* UPLO,int * N,int * NRHS,double* A,int* LDA,double* B,int* LDB,double* INFO );
           
         
"""

mkl_options = {
  "darwin_10.3_64":("%s/lib","-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "darwin_10.2_64":("%s/lib/em64t","-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core   -liomp5 -lpthread -lm"),
  "darwin_10.1_64":("%s/lib/em64t","-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "darwin_10.0_64":("%s/lib/em64t","-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.0_64" :("%s/lib/em64t","-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.1_64" :("%s/lib/em64t","-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.0_32" :("%s/lib/32","-lmkl_intel -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.1_32" :("%s/lib/32","-lmkl_intel -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.2_32" :("%s/lib/32","-lmkl_intel -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.2_64" :("%s/lib/em64t"," -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.3_64" :("%s/lib/intel64"," -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
  "linux_10.3_32" :("%s/lib/ia32"," -lmkl_intel -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lm"),
}
