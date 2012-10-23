/*
 *  errorlist.h
 *  likeli
 *
 *  Created by Karim Benabed on 06/02/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __ERRORLIST_H
#define __ERRORLIST_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


#define TXT_SZ 4192
#define WHR_SZ 2048
#define TOT_SZ TXT_SZ+WHR_SZ

typedef struct _err {
	char errWhere[WHR_SZ];
	char errText[TXT_SZ];
	int errValue;
	struct _err* next;
} error;
typedef error sm2_error;

#define forwardErr -123456789
#define undefErr   -987654321
#define noErr       0
#define placeHolderError 1020304050

#define el_alloc   -42

error* newError(int errV, char* where, char* text, error* linkTo, error* addErr);
error* newErrorVA(int errV, const char* func,char* where, char* text, error* linkTo, error* addErr,...);

void stringError(char* str, error *err);
void printError(FILE* flog,error* err);
int getErrorValue(error* err);

void purgeError(error **err); 

error* unpickleError(void* buf, int len);
int pickleError(error* err, void** buf);


error* initError(void);
void endError(error **err);

#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif



#define _DEBUGHERE_(extra,...) fprintf(stderr,"%s:("__FILE__":%d) "extra"\n",__func__,__LINE__,__VA_ARGS__);

#define someErrorVA(errV,txt,prev,next,li,...) newErrorVA(errV,__func__,"("__FILE__":"#li")",txt,prev,next,__VA_ARGS__)
#define topErrorVA(errV,txt,next,li,...)       someErrorVA(errV,txt,NULL,next,li,__VA_ARGS__)
#define oneErrorVA(errV,txt,li,...)            someErrorVA(errV,txt,NULL,NULL,li,__VA_ARGS__)
#define addErrorVA(errV,txt,prev,li,...)       someErrorVA(errV,txt,prev,NULL,li,__VA_ARGS__)

#define someError(errV,txt,prev,next,li)       someErrorVA(errV,txt,prev,next,li,"")
#define topError(errV,txt,next,li)             topErrorVA(errV,txt,next,li,"")
#define oneError(errV,txt,li)                  oneErrorVA(errV,txt,li,"")
#define addError(errV,txt,prev,li)             addErrorVA(errV,txt,prev,li,"")

// isError
#define isError(err)                         (_isError(err))

#define isErrorReturn(err,ret) {            \
  if (isError(err)) {                       \
    return ret;                             \
  }                                         \
}

#define isErrorExit(err,F) {                \
  if (isError(err)) {                       \
    printError(F,err);                      \
    exit(getErrorValue(err));               \
  }                                         \
}

// forwards
#define forwardErrorNoReturn(err,line) {    \
  if (isError(err)) {                       \
    err=topError(forwardErr,"",err,line);   \
  }                                         \
}

#define forwardError(err,line,ret) {        \
  forwardErrorNoReturn(err,line);           \
  isErrorReturn(err,ret)                    \
}

// tests 
#define testErrorVA(test,error_type,message,err,line,...) {  \
  if (test) {                                                \
    err=addErrorVA(error_type,message,err,line,__VA_ARGS__); \
  }                                                          \
}

#define testErrorRetVA(test,error_type,message,err,line,return_value,...) { \
  testErrorVA(test,error_type,message,err,line,__VA_ARGS__);                \
  isErrorReturn(err,return_value);                                          \
}

#define testErrorExitVA(test,error_type,message,err,line,...) {             \
  testErrorVA(test,error_type,message,err,line,__VA_ARGS__);                \
  isErrorExit(err,stderr);                                                  \
}

#define testError(test,error_type,message,err,line) testErrorVA(test,error_type,message,err,line,"")
#define testErrorRet(test,error_type,message,err,line,return_value) testErrorRetVA(test,error_type,message,err,line,return_value,"")
#define testErrorExit(test,error_type,message,err,line) testErrorExitVA(test,error_type,message,err,line,"")

// exit and co
#define exitOnError(err,F) isErrorExit(err,F)

#define quitOnError(err,line,F) {  \
  forwardErrorNoReturn(err,line);  \
  isErrorExit(err,F);              \
}

#define quitOnErrorStr(err,line,F,str) { \
  forwardErrorNoReturn(err,line);  \
  if (isError(err)) {                       \
    printError(F,err);                      \
    fprintf(F, "%s\n", str);                \
    exit(getErrorValue(err));               \
  }                                         \
}

#define err_base      -100
#define err_io        -1 + err_base
#define err_allocate  -2 + err_base

error* initError(void);
void endError(error **err);
int _isError(error *err);
#endif
