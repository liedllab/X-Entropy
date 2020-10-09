#include "Exceptions.h"
#include <boost/python.hpp>

void translateEmptyListError(const EmptyListError &err)
{
    PyErr_SetString(PyExc_RuntimeError, err.what());
}
