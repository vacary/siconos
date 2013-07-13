#ifdef _MSC_VER

#include <stdarg.h>
#include <stdio.h>
#include <math.h>

int c99_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
  int count = -1;

  if(size != 0)
    count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
  if(count == -1)
    count = _vscprintf(format, ap);

  return count;
}

int snprintf(char* str, size_t size, const char* format, ...)
{
  int count;
  va_list ap;

  va_start(ap, format);
  count = c99_vsnprintf(str, size, format, ap);
  va_end(ap);

  return count;
}

double rint(double x)
{
  //middle value point test
  if (ceil(x+0.5) == floor(x+0.5))
  {
    int a = (int)ceil(x);
    if (a%2 == 0)
    {return ceil(x);}
    else
    {return floor(x);}
  }

  else return floor(x+0.5);
}

#endif
