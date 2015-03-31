#ifndef _RIEMANNSOLVER_PHYSICAL_OBJECTS_H_
#define _RIEMANNSOLVER_PHYSICAL_OBJECTS_H_

#include "gas.h"
#include "tiny_stack.h"

typedef enum _Object_Type {
	shock_wave = 0,
	expansion_wave = 1,
	material_discontinuity = 2,
	trivial = 3
} Object_Type;

/* no more than 4 is needed to describe the phenomena of discontinuity problem */
#define STACK_VOLUME 4

typedef Tiny_Stack<Status, STACK_VOLUME> Status_Container;

typedef struct _Physical_Object {
	Status_Container::pointer front;
	Status_Container::pointer back;
	Object_Type type;
} Physical_Object;

#endif