#include "Exceptions.h"
#include <boost/python.hpp>

void translateEmptyListError(EmptyListError &err)
{
    PyErr_SetString(PyExc_RuntimeError, err.what());
}