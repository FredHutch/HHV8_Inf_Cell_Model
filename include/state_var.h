#ifndef STATEVAR_H
#define STATEVAR_H
class state_var
{
    private:
	unsigned long int value;

    public:
	void initialz_with (unsigned int val);
	void update_with (long int val);
	unsigned long int val ();
};
#endif
