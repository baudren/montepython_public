/*
 *  errorlist.c
 * simplified version from pmclib
 */

#ifdef __PLANCK__
#include "HL2_likely/tools/errorlist.h"
#else
#include "errorlist.h"
#endif

error* newError(int errV, char* where, char* text, error* linkTo, error* addErr) {
  error* er;
  int sts;
	
  /*fprintf(stderr,"new : %d %s %s (%p %p)\n",errV,where,text,linkTo,addErr); */
  er = (error*) malloc(sizeof(error));
  if (er == NULL) {
    fprintf(stderr,"Cannot allocate memory for error management! Aborting\n");
    exit(-1);
  }

  sts = snprintf(er->errWhere,WHR_SZ,"%s",where);
  if (sts<=0) {
    sprintf(er->errWhere,"PROBLEM IN ERROR TEXT");
  }
  if (sts==WHR_SZ) {
    er->errWhere[WHR_SZ-1]='\0';
  }
  sts = snprintf(er->errText,TXT_SZ,"%s",text);
  if (sts<=0) {
    sprintf(er->errText,"PROBLEM IN ERROR TEXT");
  }
  if (sts==TXT_SZ) {
    er->errText[TXT_SZ-1]='\0';
  }
   er->errValue = errV;
   er->next=NULL;
  if (isError(addErr)) {
    er->next = addErr;    
  } else {
    purgeError(&addErr);
  }
  if (isError(linkTo)) {
    linkTo->next = er;
    return linkTo;
  } else {
    purgeError(&linkTo);
  }
  return er;
}
	
error* newErrorVA(int errV,const char* func,char* where, char* text,error* linkTo, error* addErr,...) {
  char logTXTA[TXT_SZ],logTXT[TXT_SZ];
  va_list args;
  int sz;
  
  sprintf(logTXTA,"%s%s",func,where);
  va_start(args,addErr);
  sz=vsnprintf(logTXT, TXT_SZ,text,args);
  va_end(args);
  if (sz<=0) {
    sprintf(logTXT,"PROBLEM IN ERROR TEXT");
  }
  if (sz==TXT_SZ) {
    logTXT[TXT_SZ-1]='\0';
  }
  return newError(errV,logTXTA,logTXT,linkTo,addErr);
}

void purgeError(error **err) {
   if (*err==NULL)
     return;
   if ((*err)->next!=NULL)
     purgeError(&((*err)->next));
   free(*err);
   *err = NULL;
}


void stringError(char* str, error *err) {
   error *nerr;
   char spc[512];

   sprintf(str,"%s", "");	
   sprintf(spc,"%s", "");
	
   for(nerr=err;nerr!=NULL;nerr=nerr->next) {
    if (nerr->errValue==forwardErr) {
	sprintf(str,"%s%s%s::ForwardError\n",str,spc,nerr->errWhere);
    } else {
      sprintf(str,"%s%s%s::Error %d (%s)\n",str,spc,nerr->errWhere,nerr->errValue,nerr->errText);
    }
      sprintf(spc,"%s  ",spc);
   }
}
	
void printError(FILE* flog,error *err)
{
   char str[(TXT_SZ+WHR_SZ)*50];
   FILE* nf;
   nf=flog;
   if (flog==NULL) 
     nf=stderr;
   stringError(str,err);
   fprintf(nf,"%s",str);
}

int getErrorValue(error *err)
{
   if (err==NULL) {
      return noErr;
   }
   if (err->errValue==forwardErr) {
      if (err->next==NULL)
	return undefErr;
      return getErrorValue(err->next);
   }	
   return err->errValue;
}

int pickleError(error* err, void** buf) {
  int len;
  error* cur;
  void* bufcur;
  /* compute size */
  len=1;
  cur = err;
  while (cur->next!=NULL) {
    cur = cur->next;
    len++;
  }

  *buf = malloc(sizeof(error)*len);
  cur = err;
  bufcur = *buf;
  while(cur!=NULL) {
    memcpy((void*) bufcur,(void*) cur,sizeof(error));
    /*if (cur->next!=NULL) {
      ((error*) bufcur)->next = sizeof(error)*(len+1);
    }*/
    cur = cur->next;
    bufcur =  ((char*)bufcur)+sizeof(error); 
  }
  return len;
}

error* unpickleError(void* buf, int len) {
  int i;
  error *cur,*prev,*first;
  prev=NULL;
  first=NULL;
  for(i=0;i<len;i++) {
    cur = malloc(sizeof(error));
    memcpy(cur,((char*)buf)+sizeof(error)*i,sizeof(error));
    cur->next=NULL;
    if (prev!=NULL) {
      prev->next=cur;
    } else {
      first=cur;
    }
    prev=cur;
  }
  return first;
}

error* initError(void) {
  error *err;
  
  err = NULL;
  err=addError(placeHolderError,"PlaceHolder",NULL,0); 
  return err;
  
}

void endError(error **err) {
  purgeError(err);
}

int _isError(error *err) {
  if (err==NULL) {
    return 1==0;
  }
  if ((err)->errValue == placeHolderError) {
    return 1==0;
  }
  return 1==1;
}
