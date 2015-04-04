#ifndef _DISPLAY_H_
#define _DISPLAY_H_

#define DISPLAY_TEST 0
#define PI 3.14159265f

class Display{
public:

	void Run(int, char**);

};

void GLInit(int, char**);
void idleDisplay(void);
void reShape(int, int);
void displayVec(float, float, float, float, float, float);
void displayFunc();
void keyboardFunc(unsigned char, int, int);
void mouseButton(int, int, int, int);
void mouseMove(int, int);


#endif