/* ============================================================ *
 * io.c                                                         *
 * Martin Kilbinger, Karim Benabed                              *
 * ============================================================ */
/* simplified version from pmclib*/

#include "io.h"

/* ============================================================ *
 * I/O and memory tools.					*
 * ============================================================ */

FILE* fopen_err(const char* fname, const char *mode, error **err)
{
  FILE* fp;
  fp = fopen(fname, mode);
  testErrorRetVA(fp==NULL, err_io, "Cannot open file '%s' (mode \"%s\")", *err,__LINE__,NULL,fname, mode);
  return fp;
}

void* malloc_err(size_t sz, error **err)
{
  void* res;
  
  testErrorRetVA(sz<=0,err_allocate,"Size %ld has to be positive",*err,__LINE__,NULL,sz);
  
  res = malloc(sz);
  testErrorRetVA(res==NULL,err_allocate,"Cannot allocate %ld bytes",*err,__LINE__,NULL,sz);
  return res;
}

void* calloc_err(size_t nl, size_t sz1, error **err)
{
  void* res;
  testErrorRetVA(nl<=0,err_allocate,"Size %ld has to be positive",*err,__LINE__,NULL,nl);
  
  res = calloc(nl,sz1);
  testErrorRetVA(res==NULL,err_allocate,"Cannot allocate %d objects of %d bytes (=%d total)",
		 *err,__LINE__,NULL,nl,sz1,nl*sz1);
  return res;
}

void *realloc_err(void *ptr, size_t sz, error **err)
{
   void *res;

   res = realloc(ptr, sz);
   testErrorVA(ptr==NULL, err_allocate, "Cannot reallocate %zu bytes", *err, __LINE__, NULL, sz);
   testErrorVA(res==NULL, err_allocate, "Cannot reallocate %zu bytes", *err, __LINE__, NULL, sz);

   return NULL;
}

/* Returns a pointer with content orig, resized from sz_orig to sz_new. If *
 * free_orig=1, the original memory is freed.				   
 * orig can be NULL or sz_orig can be 0, in that case, just allocate sz_new and returns
 */
void* resize_err(void* orig, size_t sz_orig, size_t sz_new, int free_orig, error **err)
{
  void* res;
  size_t minsz,o_size;

  res = malloc_err(sz_new,err);
  forwardError(*err,__LINE__,NULL);

  if (orig == NULL || sz_orig==0) {
    return res;
  }
  
  o_size = sz_orig;
  
  
  minsz = o_size;
  if (o_size > sz_new) {
    minsz = sz_new;
  }
  memcpy(res,orig,minsz);
  if(free_orig == 1) {
    free(orig);
  }
  return res;
}

size_t read_line(FILE* fp, char* buff, int bmax, error **err) {
  size_t i;
  char cur;
  int dsc;
  
  i=0;
  dsc=0;
  while (i<bmax-1) {
     cur=(char)fgetc(fp);
    if (feof(fp) || cur=='\n') {
      buff[i]='\0';
      //fprintf(stderr,"Read |%s|\n",buff);
      return i+1;
    }
    if (dsc==1) {
      continue;
    }
    if (cur=='#' || cur=='!' || cur==';') {
      dsc=1;
      continue;
    }
    buff[i]=cur;
    i++;
  }
  *err = addErrorVA(err_io,"One line of file is longer than the max size (%d)",*err,__LINE__,bmax);
  return -1;
}


void* read_any_vector(const char* fnm, size_t n, char *prs, size_t sz, error **err) {
  size_t nn;
  char* res;
  nn=n;
  res = read_any_list(fnm, &nn, prs, sz, err);
  forwardError(*err,__LINE__,NULL);
  return res;
}


//#define _stp_ 4096
#define _stp_ 524288

void* read_any_list(const char* fnm, size_t *n, char *prs, size_t sz, error **err)
{
   void *res;
   size_t nlines;

   res = read_any_list_count(fnm, n, prs, sz, &nlines, err);
   forwardError(*err, __LINE__, NULL);

   return res;
}

/* ============================================================ *
 * Read white-spaced separated records of format prs and size   *
 * sz from the file fnm. The number of records is set in n on   *
 * output. On input, n should be set to zero, or ask Karim what *
 * happens if not. Returns record as void pointer.		*
 * The number of lines is set in nlines on output.		*
 * Comment and empty lines are ignored.			        *
 * ============================================================ */
void* read_any_list_count(const char* fnm, size_t *n, char *prs, size_t sz,
			  size_t *nlines, error **err)
{
  FILE *fp;
  size_t i,cz,j,j0,ibefore;
  char *res;
  char buffer[_stp_];
  int clip;
  
  if (*n>0) {
    clip=1;
  } else {
    clip=0;
    *n=_stp_;
  }
  
  res = (void*) malloc_err((*n)*sz,err);
  forwardError(*err,__LINE__,NULL);
  
  fp=fopen_err(fnm,"r",err);
  forwardError(*err,__LINE__,NULL);
  
  i=0;
  *nlines = 0;
  while(feof(fp)==0) {
    cz=read_line(fp,buffer,_stp_,err);
    testErrorRetVA(isError(*err),err_io,"Reading file '%s'",*err,__LINE__,NULL,fnm);
    j=0;
    //fprintf(stderr,"A %d %s\n",sz,prs);
    ibefore = i;
    while (j<cz) {
       //fprintf(stderr,"%d %c\n",j,buffer[j]);
      //look for first non blanc char
      while(j<cz && isspace(buffer[j]) && buffer[j]!='\0') {
        //fprintf(stderr,"%d %c\n",j,buffer[j]);
        j++;
      }
      if (j==cz || buffer[j]=='\0') //blank line
        break;
      // find next blanck char
      j0=j;
      while(j<cz && !isspace(buffer[j]) && buffer[j]!='\0')
        j++;
      buffer[j]='\0';
      j++;
      //fprintf(stderr,"segment |%s|\n",&(buffer[j0]));
      // read value between the two limits
      testErrorRetVA(sscanf(buffer+j0,prs,res +i*sz)!=1,err_io,"Cannot read file '%s' at index %d",
		     *err,__LINE__,NULL,fnm,i)
	//fprintf(stderr,prs,*(double*)(res +i*sz));
	//fprintf(stderr,"\n");
      i++;
      if (i==*n) {
        if (clip==1) 
          break;
        res = resize_err(res,(*n)*sz,(*n)*2*sz,1,err);
        forwardError(*err,__LINE__,NULL);
        (*n) = (*n)*2;
      }
    }

    /* NEW (MK) */
    if (i>ibefore) {
       /* A valid line was found: increase line count */
       (*nlines) ++;
    }

    if (i==*n) {
      if (clip==1) 
        break;
      res = resize_err(res,(*n)*sz,(*n)*2*sz,1,err);
      forwardError(*err,__LINE__,NULL);
      (*n) = (*n)*2;
    }
  }
  fclose(fp);
  if (clip==1)
    testErrorRetVA(i!=*n,err_io,"Not enough data in file '%s' (got %d entries, expected %d)\n",
		   *err,__LINE__,NULL,fnm,i,*n);
  
  res = resize_err(res,(*n)*sz,i*sz,1,err);
  forwardError(*err,__LINE__,NULL);
  *n=i;
  return res;
}
#undef _stp_

double* read_double_vector(const char* fnm, size_t n, error **err) {
  double *res;
  res = read_any_vector(fnm, n, "%lg", sizeof(double), err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

double* read_double_list(const char* fnm, size_t *n, error **err) {
  double *res;
  res = read_any_list(fnm, n, "%lg", sizeof(double), err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

float* read_float_vector(const char* fnm, size_t n, error **err) {
  float *res;
  res = read_any_vector(fnm, n, "%g", sizeof(float), err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

float* read_float_list(const char* fnm,size_t *n, error **err) {
  float *res;
  res = read_any_list(fnm, n, "%g",sizeof(float),err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

int* read_int_vector(const char* fnm,size_t n, error **err) {
  int *res;
  res = read_any_vector(fnm, n, "%d",sizeof(int),err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

int* read_int_list(const char* fnm, size_t *n, error **err) {
  int *res;
  res = read_any_list(fnm, n, "%d",sizeof(int),err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

long* read_long_vector(const char* fnm, size_t n, error **err) {
  long *res;
  res = read_any_vector(fnm, n, "%ld",sizeof(long),err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

long* read_long_list(const char* fnm, size_t *n, error **err) {
  long *res;
  res = read_any_list(fnm, n, "%ld",sizeof(long),err);
  forwardError(*err,__LINE__,NULL);
  return res;
}

void* read_bin_vector(const char* fnm, size_t n, error **err) {
  void *res;
  FILE* ff;
  size_t rr;
  
  ff = fopen_err(fnm, "r", err);
  forwardError(*err,__LINE__,NULL);
  res = malloc_err(n,err);
  forwardError(*err,__LINE__,NULL);
  rr=fread(res, 1, n, ff);
  testErrorRetVA(rr!=n,io_file,"File '%s' does not contain enough data (Expected %d bytes, got %d)",*err,__LINE__,NULL,fnm,n,rr);
  fclose(ff);
  return res;
}

void* read_bin_list(const char* fnm, size_t *n, error **err) {
  void *res;
  FILE* ff;
  size_t rr;
  struct stat stst;
  
  testErrorRetVA(stat(fnm,&stst)!=0,io_file,"File '%s' cannot be stated",*err,__LINE__,NULL,fnm);
  *n = stst.st_size;
  
  ff = fopen_err(fnm, "r", err);
  forwardError(*err,__LINE__,NULL);
  res = malloc_err(*n,err);
  forwardError(*err,__LINE__,NULL);
  rr=fread(res, 1, *n, ff);
  testErrorRetVA(rr!=*n,io_file,
		 "File '%s' does not contain enough data (Expected %d bytes, got %d)",
		 *err,__LINE__,NULL,fnm,*n,rr);
  fclose(ff);
  return res;
}

void write_bin_vector(void* buf, const char* fnm, size_t n, error **err) {
  FILE* ff;
  size_t rr;
  
  ff = fopen_err(fnm, "w", err);
  forwardError(*err,__LINE__,);
  rr=fwrite(buf, 1, n, ff);    
  
  testErrorRetVA(rr!=n,io_file,
		 "Cannot write data in file '%s' (Expected %d bytes, got %d)",
		 *err,__LINE__,,fnm,n,rr);
  fclose(ff);
  fprintf(stderr,"Saved in %s\n", fnm);
}

time_t start_time(FILE *FOUT)
{
   time_t t_start;
   time(&t_start);
   fprintf(FOUT, "Started at %s\n", ctime(&t_start));
   fflush(FOUT);
   return t_start;
}

void end_time(time_t t_start, FILE *FOUT)
{
   time_t t_end;
   double diff;

   time(&t_end); 
   fprintf(FOUT, "Ended at %s", ctime(&t_end));
   diff = difftime(t_end, t_start);
   fprintf(FOUT, "Computation time %.0fs (= %dd %dh %dm %ds)\n",
           diff, (int)diff/86400, ((int)diff%86400)/3600,
           ((int)diff%3600)/60, ((int)diff%60));
}

