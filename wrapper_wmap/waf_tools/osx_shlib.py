def configure(ctx):
  import sys
  ctx.env.shsuffix = "so"
  if sys.platform.lower()=="darwin":
    ctx.env.shsuffix = "dylib"