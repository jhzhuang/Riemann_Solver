#ifndef _RIEMANNSOLVER_TINY_STACK_H_
#define _RIEMANNSOLVER_TINY_STACK_H_

template<typename T, int Volume>
class Tiny_Stack {
public:
	typedef int pointer;
public:
	Tiny_Stack();
	~Tiny_Stack() {}
public:
	inline void push(T &val);
	inline T& pop(void);
	inline T& at(pointer ptr);
	inline pointer get_hwm();
private:
	T stack[Volume];
	int hwm;    //high water mark pointer
};

template<typename T, int Volume>
Tiny_Stack<T, Volume>::Tiny_Stack(): hwm(0) {}

template<typename T, int Volume>
void Tiny_Stack<T, Volume>::push(T &val) {
	if (hwm < Volume)
	{
		stack[hwm] = val;
		hwm = hwm + 1;
	}
}

template<typename T, int Volume>
T& Tiny_Stack<T, Volume>::pop() {
	if (hwm > -1)
	{
		hwm = hwm - 1;
		return stack[hwm];
	}
}

template<typename T, int Volume>
T& Tiny_Stack<T, Volume>::at(pointer ptr) {
	return stack[ptr];
}

template<typename T, int Volume>
typename Tiny_Stack<T, Volume>::pointer Tiny_Stack<T, Volume>::get_hwm() {
	return hwm;
}
#endif