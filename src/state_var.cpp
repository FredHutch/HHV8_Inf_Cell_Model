#include "state_var.h"

void state_var::initialz_with (unsigned int val)
{
 value = val;
} 
void state_var::update_with (long int val)
{
    if (val > 0 || value > -val)
      value += val;
    else
      value = 0;
}
unsigned long int state_var::val ()
{
  return (value);
}
