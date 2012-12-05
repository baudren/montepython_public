class LockException(Exception):
  # Error codes:
  LOCK_FAILED = 1

def lock(file, flags):
  import fcntl
  try:
    fcntl.flock(file.fileno(), flags)
  except IOError, exc_value:
    # The exception code varies on different systems so we'll catch
    # every IO error
    raise LockException(*exc_value)

def unlock(file):
  import fcntl
  fcntl.flock(file.fileno(), fcntl.LOCK_UN)

