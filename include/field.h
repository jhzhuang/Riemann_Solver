#ifndef _RIEMANNSOLVER_FIELD_H_
#define _RIEMANNSOLVER_FIELD_H_

#include "physical_objects.h"

typedef struct _Field {
	Physical_Object left;
//	Physical_Object central;
	Physical_Object right;
} Field;

#endif