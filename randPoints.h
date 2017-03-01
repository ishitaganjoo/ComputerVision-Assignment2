#ifndef __RANDPOINTS_H__
#define __RANDPOINTS_H__


class randPoints{


public:
	
    float _x;
    float _y;
    float _xP;
    float _yP;

randPoints() {
  }

  randPoints(float x, float y, float xP, float yP){
    _x = x;
    _y = y;
    _xP  = xP;
    _yP = yP;
  }

};

#endif
